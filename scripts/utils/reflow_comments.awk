#!/usr/bin/awk -f
# Move over-long trailing '# ...' comments onto their own line(s) ABOVE the code
# they annotate, wrapped to MAXLEN columns. This is the repo's comment
# convention: a line whose trailing '# what + why' would push it past 100
# columns gets the comment hoisted above instead, so nothing is lost to a
# terminal's right edge and no explanation has to be shortened.
#
# Reads one script on stdin (or as $1) and writes the reflowed version to
# stdout -- it never edits in place, so the caller decides what to keep:
#     awk -f scripts/utils/reflow_comments.awk foo.sh > foo.new && mv foo.new foo.sh
# or over a whole tree:
#     for f in scripts/**/*.sh; do
#         awk -f scripts/utils/reflow_comments.awk "$f" > "$f.new" \
#             && bash -n "$f.new" && mv "$f.new" "$f"
#     done
#
# What it will NOT touch, on purpose:
#   * comments inside multi-line single-quoted blocks (the embedded awk
#     programs in kraken2_filter_reads*.sh, extract_u3_by_hxb2_anchor.sh,
#     run_hivseqinr.sh, setup_jaspar.sh, ...) -- quote state is carried across
#     lines and such lines are emitted verbatim. Hoisting a comment out of
#     another language's syntax is not safe, so those stay trailing.
#   * heredoc bodies -- emitted verbatim.
#   * lines already at or under MAXLEN -- their trailing comment stays put,
#     which is why the result is intentionally a mix of both styles.
#
# It is also careful that:
#   * a '#' only starts a comment when unquoted AND at the start of a word, so
#     "$#", "${A[@]}", "^#\|x" and 'a # b' are never mistaken for comments
#   * a comment on a backslash-continuation line is hoisted above the FIRST
#     physical line of that command, never inserted mid-continuation
#   * running it twice changes nothing (idempotent)
#
# ALWAYS verify a bulk run two ways, not one. 'bash -n' alone is not enough --
# it still passes if a comment was wrongly lifted out of a quoted string. Also
# diff the comment-stripped "code stream" before/after and require zero
# difference. See feedback_comment_style_above_line in the project memory.

BEGIN { MAXLEN = 100; n = 0; in_hd = 0; cont = 0; cmd_start = 0; PSQ = 0; PDQ = 0 }

# Walk a line, carrying PSQ/PDQ (single/double quote state) across calls.
# Returns the position of the comment '#', or 0 when the line has none.
# Only meaningful for lines that began unquoted.
function scan_line(s,    i, c, p, len) {
    len = length(s)
    for (i = 1; i <= len; i++) {
        c = substr(s, i, 1)
        if (PSQ) { if (c == "'") PSQ = 0; continue }   # inside '...': only ' ends it
        if (PDQ) {                                    # inside "...": \ escapes
            if (c == "\\") { i++; continue }
            if (c == "\"") PDQ = 0
            continue
        }
        if (c == "\\") { i++; continue }               # escaped char, skip it
        if (c == "'") { PSQ = 1; continue }
        if (c == "\"") { PDQ = 1; continue }
        if (c == "#") {
            if (i == 1) return i                      # whole line is a comment
            p = substr(s, i - 1, 1)
            if (p == " " || p == "\t") return i        # word start => rest is a comment
        }
    }
    return 0
}

# Delimiter of a heredoc opened on this code line, else "". Sets HD_TABS for <<-.
function heredoc_delim(code,    i, len, c, rest, s2, d2) {
    HD_TABS = 0
    len = length(code); s2 = 0; d2 = 0
    for (i = 1; i <= len - 1; i++) {
        c = substr(code, i, 1)
        if (s2) { if (c == "'") s2 = 0; continue }
        if (d2) {
            if (c == "\\") { i++; continue }
            if (c == "\"") d2 = 0
            continue
        }
        if (c == "\\") { i++; continue }
        if (c == "'") { s2 = 1; continue }
        if (c == "\"") { d2 = 1; continue }
        if (c != "<") continue
        if (substr(code, i + 1, 1) != "<") continue
        if (substr(code, i + 2, 1) == "<") { i += 2; continue }   # <<< is a herestring
        rest = substr(code, i + 2)
        if (substr(rest, 1, 1) == "-") { HD_TABS = 1; rest = substr(rest, 2) }
        sub(/^[ \t]+/, "", rest)
        if (match(rest, /^"[^"]+"/) || match(rest, /^'[^']+'/)) return substr(rest, 2, RLENGTH - 2)
        if (match(rest, /^[A-Za-z_][A-Za-z_0-9]*/)) return substr(rest, 1, RLENGTH)
        i += 1
    }
    return ""
}

# Split comment text into WRAPPED[1..WN], each fitting MAXLEN once indent + "# " is added.
function wrap(text, indent,    width, words, m, i, cur, w) {
    width = MAXLEN - length(indent) - 2
    if (width < 24) width = 24                        # deeply indented code: allow overflow
    WN = 0; cur = ""
    m = split(text, words, /[ \t]+/)
    for (i = 1; i <= m; i++) {
        w = words[i]
        if (cur == "") cur = w
        else if (length(cur) + 1 + length(w) <= width) cur = cur " " w
        else { WRAPPED[++WN] = cur; cur = w }         # this word would overflow, start a line
    }
    if (cur != "") WRAPPED[++WN] = cur
    return WN
}

# Make room for k lines at position pos in out[].
function shift_down(pos, k,    i) {
    for (i = n; i >= pos; i--) out[i + k] = out[i]
}

{
    line = $0

    if (in_hd) {                                      # heredoc body: never rewritten
        t = line
        if (HD_TABS) sub(/^\t+/, "", t)
        out[++n] = line
        if (t == hd_delim) in_hd = 0                  # terminator reached
        next
    }

    started_quoted = (PSQ || PDQ)                     # did this line open inside a string?
    cp = scan_line(line)
    ends_quoted = (PSQ || PDQ)

    if (started_quoted) {                             # embedded program/string: verbatim
        out[++n] = line
        cont = (!ends_quoted && line ~ /\\$/) ? 1 : 0
        next
    }

    if (cp > 0) { code = substr(line, 1, cp - 1); ctext = substr(line, cp + 1) }
    else        { code = line; ctext = "" }

    codetrim = code
    sub(/[ \t]+$/, "", codetrim)

    if (!cont) cmd_start = n + 1                      # first physical line of this command

    if (cp > 0 && length(line) > MAXLEN && codetrim !~ /^[ \t]*$/) {
        target = (cmd_start <= n) ? cmd_start : 0     # 0 = this line is its own command start
        ind = (target > 0) ? out[target] : codetrim   # indent comes from the target line
        match(ind, /^[ \t]*/); indent = substr(ind, 1, RLENGTH)

        gsub(/^[ \t]+|[ \t]+$/, "", ctext)
        wrap(ctext, indent)

        if (target > 0) {                             # hoist above the command's first line
            shift_down(target, WN)
            for (k = 1; k <= WN; k++) out[target + k - 1] = indent "# " WRAPPED[k]
            n += WN
        } else {
            for (k = 1; k <= WN; k++) out[++n] = indent "# " WRAPPED[k]
        }
        out[++n] = codetrim                           # the code, comment now stripped
    } else {
        keep = line
        sub(/[ \t]+$/, "", keep)                      # drop trailing whitespace while here
        out[++n] = keep
    }

    # Continuation + heredoc state, judged on the code half only, and only when
    # the line ended outside any quote (otherwise the command continues inside a
    # string and the next line is passed through verbatim anyway).
    if (ends_quoted) {
        cont = 0
    } else {
        cont = (codetrim ~ /\\$/) && (codetrim !~ /\\\\$/)
        d = heredoc_delim(codetrim)
        if (d != "") { in_hd = 1; hd_delim = d }
    }
}

END { for (i = 1; i <= n; i++) print out[i] }
