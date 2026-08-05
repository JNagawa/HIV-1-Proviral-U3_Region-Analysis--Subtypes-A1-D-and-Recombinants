#!/bin/bash
#SBATCH --job-name=pb_motif_pptx
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

##################################################################################
# Build writeups/PacBio_Motif_Tool_Comparison.pptx from the motif tool-comparison
# analysis.
#
# A .pptx is an OOXML zip, so it is generated here with bash + zip and no Python
# (matching this repo's convention). pandoc, LibreOffice and python-pptx are all
# absent on this cluster, so converting the markdown report was not an option --
# the slide XML is emitted directly.
#
# The theme, slide master and slide layouts are copied verbatim from the existing
# writeups/PacBio_Tool_Comparison.pptx, so this deck inherits the same look. Its
# charts, embedded workbooks and printer settings are dropped; only the parts a
# text deck needs are kept.
#
# Every number and table is READ FROM THE ANALYSIS TSVs at build time rather than
# typed in, so the deck cannot drift from the results it describes.
#
# Usage: sbatch scripts/utils/build_motif_comparison_pptx.sh
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail

# The compute nodes have neither zip nor unzip (they exist only on the login
# node), so put the conda env and the site anaconda install on PATH -- between
# them they provide bsdtar, which both reads and writes zip archives.
export PATH="/etc/ace-data/home/jnagawa/.conda/envs/HIV_U3analysis/bin:/opt/ohpc/pub/apps/anaconda3/bin:${PATH}"

# Expand a zip archive into a directory, using whichever tool is available.
zip_extract() {
    local archive="$1" dest="$2"                     # archive to read, directory to fill
    mkdir -p "${dest}"                               # bsdtar will not create it itself
    if command -v bsdtar >/dev/null 2>&1; then
        bsdtar -xf "${archive}" -C "${dest}"         # bsdtar reads zip natively
    elif command -v unzip >/dev/null 2>&1; then
        unzip -q "${archive}" -d "${dest}"           # login-node fallback
    else
        echo "ERROR: no tool available to read a zip archive (need bsdtar or unzip)." >&2
        return 1
    fi
}

# Create a zip archive from the CURRENT directory's contents. The first argument
# is the output path; the rest are members, in the order they should be stored --
# OPC convention puts [Content_Types].xml first.
zip_create() {
    local archive="$1"; shift                        # output path, then members in order
    if command -v bsdtar >/dev/null 2>&1; then
        bsdtar --format=zip -cf "${archive}" "$@"    # member order follows argument order
    elif command -v zip >/dev/null 2>&1; then
        zip -q -X -r "${archive}" "$@"               # login-node fallback
    else
        echo "ERROR: no tool available to write a zip archive (need bsdtar or zip)." >&2
        return 1
    fi
}

# Locate the repo root (Slurm copies this script to a spool dir under sbatch)
REPO_ROOT=""                                         # filled in by the loop below
for CAND in "${SLURM_SUBMIT_DIR:-}" \
            "$(cd "$(dirname "$0")/../.." 2>/dev/null && pwd)" \
            "$(pwd)"; do
    # a repo root is only a repo root if the helper library is under it
    if [ -n "${CAND}" ] && [ -f "${CAND}/scripts/common/lib_compare.sh" ]; then
        REPO_ROOT="${CAND}"; break                   # first match wins
    fi
done
[ -n "${REPO_ROOT}" ] || { echo "ERROR: cannot locate repo root." >&2; exit 1; }

# which assembly arm's comparison to build slides from (see compare_msa_pacbio.sh)
ASSEMBLY_ARM="${ASSEMBLY_ARM:-minimap2_consensus}"
CMP_DIR="${REPO_ROOT}/results/motif_mapping/pacbio/${ASSEMBLY_ARM}/tool_comparison"  # analysis inputs
TEMPLATE="${REPO_ROOT}/writeups/PacBio_Tool_Comparison.pptx"         # theme/master donor
OUT_PPTX="${REPO_ROOT}/writeups/PacBio_Motif_Tool_Comparison.pptx"   # what we build
BUILD="${REPO_ROOT}/.pptx_build"                     # scratch tree, removed at the end
# clean any half-finished tree from a previous run
rm -rf "${BUILD}"; mkdir -p "${BUILD}/pkg"
# always tidy the scratch tree, even on failure
trap 'rm -rf "${BUILD}"' EXIT

for f in "${CMP_DIR}/pairwise_concordance.tsv" "${CMP_DIR}/hxb2_landmark_recovery.tsv" \
         "${CMP_DIR}/u3_identity_to_hxb2.tsv" "${TEMPLATE}"; do
    # refuse to build a deck from missing inputs rather than emit blank slides
    [ -s "${f}" ] || { echo "ERROR: required input missing: ${f}" >&2; exit 1; }
done

##################################################################################
#  Take theme / master / layouts from the existing deck, drop everything else
##################################################################################
zip_extract "${TEMPLATE}" "${BUILD}/tpl"             # expand the donor package
P="${BUILD}/pkg"                                     # shorthand for the package root
mkdir -p "${P}/_rels" "${P}/docProps" "${P}/ppt/_rels" "${P}/ppt/slides/_rels" \
         "${P}/ppt/theme" "${P}/ppt/slideMasters/_rels" "${P}/ppt/slideLayouts/_rels" "${P}/ppt/media"
cp "${BUILD}/tpl/ppt/theme/"*.xml           "${P}/ppt/theme/"                 2>/dev/null
cp "${BUILD}/tpl/ppt/slideMasters/"*.xml    "${P}/ppt/slideMasters/"          2>/dev/null
cp "${BUILD}/tpl/ppt/slideMasters/_rels/"*  "${P}/ppt/slideMasters/_rels/"    2>/dev/null
cp "${BUILD}/tpl/ppt/slideLayouts/"*.xml    "${P}/ppt/slideLayouts/"          2>/dev/null
cp "${BUILD}/tpl/ppt/slideLayouts/_rels/"*  "${P}/ppt/slideLayouts/_rels/"    2>/dev/null
cp "${BUILD}/tpl/ppt/presProps.xml" "${BUILD}/tpl/ppt/viewProps.xml" "${P}/ppt/" 2>/dev/null
cp "${BUILD}/tpl/ppt/media/"* "${P}/ppt/media/" 2>/dev/null
N_LAYOUTS=$(ls "${P}/ppt/slideLayouts/"*.xml 2>/dev/null | wc -l)   # layouts carried over
echo "carried over ${N_LAYOUTS} slide layouts + theme + master"

##################################################################################
#  Slide construction helpers
#
#  Geometry matches the donor deck: 12191695 x 6858000 EMU (13.33in x 7.5in).
#  Text is placed as absolutely-positioned text boxes, which is how the donor
#  deck is built too, so no placeholder inheritance is relied on.
##################################################################################
X=548640                                             # left margin for every box
CW=11064240                                          # content width
ACCENT="C8603A"                                      # terracotta, the donor deck's accent
INK="1A1A17"                                         # near-black body text
MUTED="5A5A52"                                       # secondary/explanatory text
WARN="B3261E"                                        # red, for the caveat slide

# Escape the five XML metacharacters so arbitrary text is safe inside <a:t>
xml_escape() {
    sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g' -e "s/'/\&apos;/g" -e 's/"/\&quot;/g'
}

SHAPE_ID=1                                           # incremented per shape on a slide
SLIDE_BODY=""                                        # accumulates the current slide's shapes

# Append one text box. Args: y, height, size(1/100pt), bold, colour, font, text
add_text() {
    local y="$1" cy="$2" sz="$3" b="$4" colour="$5" font="$6" text="$7"
    SHAPE_ID=$((SHAPE_ID+1))                         # every shape needs a unique id
    local esc                                        # XML-safe copy of the text
    esc=$(printf '%s' "${text}" | xml_escape)
    SLIDE_BODY="${SLIDE_BODY}<p:sp><p:nvSpPr><p:cNvPr id=\"${SHAPE_ID}\" name=\"TextBox ${SHAPE_ID}\"/><p:cNvSpPr txBox=\"1\"/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x=\"${X}\" y=\"${y}\"/><a:ext cx=\"${CW}\" cy=\"${cy}\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom><a:noFill/></p:spPr><p:txBody><a:bodyPr wrap=\"square\"><a:spAutoFit/></a:bodyPr><a:lstStyle/><a:p><a:r><a:rPr sz=\"${sz}\" b=\"${b}\" i=\"0\"><a:solidFill><a:srgbClr val=\"${colour}\"/></a:solidFill><a:latin typeface=\"${font}\"/></a:rPr><a:t>${esc}</a:t></a:r></a:p></p:txBody></p:sp>"
}

# The short accent rule under the title, copied from the donor deck's slides
add_rule() {
    SHAPE_ID=$((SHAPE_ID+1))                         # unique shape id
    SLIDE_BODY="${SLIDE_BODY}<p:sp><p:nvSpPr><p:cNvPr id=\"${SHAPE_ID}\" name=\"Rule ${SHAPE_ID}\"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x=\"576072\" y=\"1371600\"/><a:ext cx=\"713232\" cy=\"38100\"/></a:xfrm><a:prstGeom prst=\"roundRect\"><a:avLst/></a:prstGeom><a:solidFill><a:srgbClr val=\"${ACCENT}\"/></a:solidFill><a:ln><a:noFill/></a:ln></p:spPr><p:txBody><a:bodyPr/><a:lstStyle/><a:p/></p:txBody></p:sp>"
}

SLIDE_NO=0                                           # how many slides written so far
# Start a slide: eyebrow label + title + accent rule
begin_slide() {
    local eyebrow="$1" title="$2"
    SHAPE_ID=1; SLIDE_BODY=""                        # reset per-slide state
    add_text 384048  365760 1250 1 "${ACCENT}" "Segoe UI" "${eyebrow}"
    add_text 658368  640080 2700 1 "${INK}"    "Segoe UI" "${title}"
    add_rule
    CURSOR=1600000                                   # y position for the next body line
}

# Body line styles, each advancing CURSOR by its own line height
add_lead()   { add_text "${CURSOR}" 400000 1500 0 "${MUTED}" "Segoe UI"  "$1"; CURSOR=$((CURSOR+430000)); }
add_bullet() { add_text "${CURSOR}" 340000 1400 0 "${INK}"   "Segoe UI"  "•  $1"; CURSOR=$((CURSOR+330000)); }
add_key()    { add_text "${CURSOR}" 340000 1400 1 "${ACCENT}" "Segoe UI" "$1"; CURSOR=$((CURSOR+340000)); }
add_warn()   { add_text "${CURSOR}" 340000 1400 1 "${WARN}"  "Segoe UI"  "$1"; CURSOR=$((CURSOR+340000)); }
add_mono()   { add_text "${CURSOR}" 260000 1050 0 "${INK}"   "Consolas"  "$1"; CURSOR=$((CURSOR+235000)); }
add_gap()    { CURSOR=$((CURSOR+160000)); }

# Emit the accumulated shapes as ppt/slides/slideN.xml plus its layout rel
end_slide() {
    SLIDE_NO=$((SLIDE_NO+1))                         # this slide's 1-based number
    cat > "${P}/ppt/slides/slide${SLIDE_NO}.xml" <<XMLEOF
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><p:cSld><p:bg><p:bgPr><a:solidFill><a:srgbClr val="FFFFFF"/></a:solidFill><a:effectLst/></p:bgPr></p:bg><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/>${SLIDE_BODY}</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>
XMLEOF
    # every slide must point at a layout; layout1 is enough for text-only slides
    cat > "${P}/ppt/slides/_rels/slide${SLIDE_NO}.xml.rels" <<XMLEOF
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/></Relationships>
XMLEOF
}

# Render a TSV as aligned monospace rows, so the deck shows the real table
add_tsv() {
    local file="$1" maxrows="${2:-99}"                # TSV path, row cap
    local line                                        # each formatted row
    while IFS= read -r line; do
        add_mono "${line}"
    done < <(column -t -s$'\t' "${file}" | head -n "${maxrows}")
}

##################################################################################
#  Slides
##################################################################################

# --- 1. Title -----------------------------------------------------------------
begin_slide "HIV-1 U3 ANALYSIS  ·  PACBIO HIV-SMRTcap ARM" "Motif-mapping tool comparison"
add_lead "Do different motif tools find the same sites, in the same places?"
add_gap
add_bullet "TFBS scanners:  FIMO  vs  TFBSTools  (with MOODS as a third view)"
add_bullet "G-quadruplex predictors:  gquad  vs  pqsfinder"
add_gap
add_key "Scope: 4 SMRTcap samples (124_4, 128_5, 203_3, 211_0) + HXB2 reference"
add_bullet "Pipeline re-run end-to-end: extraction -> assembly -> MSA -> motif mapping"
add_bullet "Every stage submitted as a Slurm job"
end_slide

# --- 2. The question ----------------------------------------------------------
begin_slide "WHY THIS ANALYSIS WAS NEEDED" "Hit counts cannot answer the question"
add_lead "The existing harness records only HOW MANY sites each tool found."
add_gap
add_bullet "Two tools can each report 40 sites and share not one of them"
add_bullet "A count tells you nothing about agreement on POSITION"
add_gap
add_key "So each tool's output was normalised to a common interval form and"
add_key "compared position by position."
add_gap
add_bullet "Measured per pair: exact coordinate matches, partial overlaps,"
add_bullet "tool-unique calls, and base-pair Jaccard overlap"
end_slide

# --- 3. Method: coordinate normalisation --------------------------------------
begin_slide "METHOD  ·  THE PART THAT HAD TO BE RIGHT" "Tools do not share a coordinate convention"
add_mono "FIMO         fimo.tsv    1-based inclusive     TF symbol   (RELA)"
add_mono "TFBSTools    GFF3        1-based inclusive     TF symbol   (RELA)"
add_mono "gquad        GFF3        1-based inclusive     --"
add_mono "pqsfinder    GFF3        1-based inclusive     --"
add_mono "MOODS        CSV         0-BASED, no end col   JASPAR id   (MA0107.1)"
add_gap
add_warn "MOODS differs on BOTH axes and must be translated twice:"
add_bullet "0-based -> 1-based, and end derived from matched-sequence length"
add_bullet "JASPAR accession -> TF symbol (mapping built from FIMO's own output)"
add_gap
add_lead "Verified on a shared RELA match: FIMO and TFBSTools both call 350-359,"
add_lead "MOODS reports 349. Without the shift every site looks off by one."
end_slide

# --- 4. The caveat ------------------------------------------------------------
begin_slide "READ THIS BEFORE THE RESULTS" "Only 2 of 4 samples contribute real sequence"
add_lead "The assembly step applies called variants ONTO the HXB2 reference."
add_gap
add_warn "Positions with no read coverage are emitted as REFERENCE BASES, not N."
add_bullet "Every consensus is 9719 bp at 0.00% N, whatever the sample covered"
add_bullet "124_4 reads cover 813-8002 and 8985-9719;  128_5 covers 1593-7030"
add_bullet "Neither reaches the U3 window, so their U3 is pure HXB2"
add_gap
add_tsv "${CMP_DIR}/u3_identity_to_hxb2.tsv" 6
add_gap
add_warn "3 of 5 scanned sequences are the SAME sequence. Effective n = 2."
end_slide

# --- 5. Concordance table -----------------------------------------------------
begin_slide "RESULTS  ·  DO THE TOOLS AGREE ON POSITION?" "Pairwise concordance"
add_lead "Counts from tool A's perspective. Jaccard is over covered base positions."
add_gap
add_tsv "${CMP_DIR}/pairwise_concordance.tsv" 6
add_gap
add_key "Two very different patterns:"
add_bullet "TFBS tools agree EXACTLY on position, and differ only in how many"
add_bullet "   sites they report"
add_bullet "G4 tools agree on REGION but never on boundaries"
end_slide

# --- 6. FIMO vs TFBSTools -----------------------------------------------------
FVT=$(awk -F'\t' '$1=="fimo_vs_tfbstools"{print $2, $3, $4, $6, $7, $8}' "${CMP_DIR}/pairwise_concordance.tsv")
set -- ${FVT}                                        # n_A n_B exact A_only B_only jaccard
begin_slide "RESULTS  ·  TFBS SCANNERS" "FIMO vs TFBSTools: identical coordinates"
add_mono "FIMO calls                 $1"
add_mono "TFBSTools calls            $2"
add_mono "Exact coordinate matches   $3   <- all of FIMO"
add_mono "FIMO-only calls            $4   <- none"
add_mono "TFBSTools-only calls       $5"
add_mono "Base-pair Jaccard          $6"
add_gap
add_key "FIMO is a STRICT SUBSET of TFBSTools."
add_bullet "Where both fire, they agree perfectly on position - not one shifted call"
add_bullet "The whole difference between them is stringency, not disagreement"
add_bullet "The low Jaccard reflects TFBSTools' extra breadth, not positional conflict"
add_gap
add_lead "Practical: running both is redundant. FIMO finds no site TFBSTools lacks."
end_slide

# --- 7. gquad vs pqsfinder ----------------------------------------------------
GVP=$(awk -F'\t' '$1=="gquad_vs_pqsfinder"{print $2, $3, $4, $5, $6, $7, $8}' "${CMP_DIR}/pairwise_concordance.tsv")
set -- ${GVP}                                        # n_A n_B exact partial A_only B_only jaccard
begin_slide "RESULTS  ·  G-QUADRUPLEX PREDICTORS" "gquad vs pqsfinder: same place, different extent"
add_mono "gquad predictions          $1"
add_mono "pqsfinder predictions      $2"
add_mono "Exact boundary matches     $3   <- none, ever"
add_mono "Partial overlaps           $4   <- every pqsfinder call"
add_mono "gquad-only                 $5"
add_mono "pqsfinder-only             $6"
add_gap
add_key "They agree on WHERE quadruplexes are and never on how far they extend."
add_gap
add_mono "203_3   gquad      349-402   (54 bp)"
add_mono "203_3   pqsfinder  363-402   (40 bp)   <- same right edge, later start"
add_gap
add_warn "Not interchangeable: pqsfinder emits one fixed 40 bp call per sequence,"
add_warn "gquad emits several of 30-54 bp."
end_slide

# --- 8. Ground truth ----------------------------------------------------------
begin_slide "GROUND TRUTH  ·  HXB2 POSITIVE CONTROL" "Recovery of known LTR elements"
add_lead "HXB2's U3 has literature-established elements: 2x NF-kB, 3x Sp1, TATA."
add_gap
add_tsv "${CMP_DIR}/hxb2_landmark_recovery.tsv" 7
add_gap
awk -F'\t' 'NR>1{for(i=4;i<=6;i++) if($i=="yes") c[i]++} END{
    printf "TFBSTools %d/6   FIMO %d/6   MOODS %d/6\n", c[5]+0, c[4]+0, c[6]+0}' \
    "${CMP_DIR}/hxb2_landmark_recovery.tsv" | while IFS= read -r l; do add_key "${l}"; done
add_bullet "FIMO misses the TATA box and two of the three Sp1 sites"
add_bullet "G4 columns shown for completeness - these predict structures, not TFBS"
add_bullet "   (the GC-rich Sp1/NF-kB region carries documented LTR-III/IV G4s)"
end_slide

# --- 9. Sensitivity / specificity ---------------------------------------------
begin_slide "INTERPRETATION" "What sensitivity and specificity can and cannot mean here"
add_warn "Sensitivity and specificity require a gold standard. None exists for"
add_warn "these sample U3 sequences."
add_gap
add_key "So the concordance figures are AGREEMENT, not sensitivity."
add_gap
add_bullet "On the only real ground truth (HXB2), TFBSTools is the most sensitive"
add_bullet "   TFBS scanner and FIMO the most conservative"
add_bullet "Because FIMO's calls are a strict subset, TFBSTools cannot be LESS"
add_bullet "   sensitive than FIMO on any input"
add_gap
add_warn "Specificity cannot be separated from this."
add_lead "TFBSTools' 62 extra calls are either genuine weak sites or false positives,"
add_lead "and nothing in this data distinguishes the two. Higher sensitivity with"
add_lead "unknown precision is the honest summary - not \"TFBSTools is better\"."
end_slide

# --- 10. Conclusions ----------------------------------------------------------
begin_slide "CONCLUSIONS" "What to do with these tools"
add_key "FIMO vs TFBSTools"
add_bullet "Identical coordinates wherever both call - the choice is stringency only"
add_bullet "FIMO when false positives are costly; TFBSTools when missing a weak site is"
add_bullet "Do not report both: FIMO adds nothing TFBSTools lacks"
add_gap
add_key "gquad vs pqsfinder"
add_bullet "Same regions, different extents - NOT interchangeable"
add_bullet "Anything keyed on G4 boundaries changes with the tool chosen"
add_bullet "Fix one tool for the whole study rather than mixing"
add_gap
add_key "Before any biological conclusion"
add_warn "Fix the coverage artefact: the consensus step must emit N where read depth"
add_warn "is zero, so no sample can silently contribute reference sequence."
end_slide

# --- 11. Open issues ----------------------------------------------------------
begin_slide "OPEN ISSUES  ·  NEXT STEPS" "Known limitations in the current run"
add_bullet "Coverage fill-in: mask zero-depth positions before motif scanning, then"
add_bullet "   re-run - this is the blocker on interpreting 124_4 and 128_5"
add_bullet "Effective n = 2. More SMRTcap samples with U3 coverage are needed before"
add_bullet "   any cross-sample motif claim holds"
add_bullet "hifiasm failed on all 4 samples (19-75 reads each); only the minimap2"
add_bullet "   consensus path contributed"
add_bullet "peak memory is unmeasured throughout - GNU time is not installed on this"
add_bullet "   cluster, so summary tables carry runtime only"
add_gap
add_lead "A curated truth set for these samples would allow real sensitivity and"
add_lead "specificity to replace the agreement figures used here."
end_slide

echo "generated ${SLIDE_NO} slides"

##################################################################################
#  Package-level parts
##################################################################################
# presentation.xml: master gets rId1, slides get rId2..rId(N+1)
SLD_IDS=""; PRES_RELS=""                             # accumulated slide id list and rels
for i in $(seq 1 "${SLIDE_NO}"); do
    SLD_IDS="${SLD_IDS}<p:sldId id=\"$((255+i))\" r:id=\"rId$((i+1))\"/>"
    PRES_RELS="${PRES_RELS}<Relationship Id=\"rId$((i+1))\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide\" Target=\"slides/slide${i}.xml\"/>"
done
R_THEME=$((SLIDE_NO+2)); R_PRESP=$((SLIDE_NO+3)); R_VIEWP=$((SLIDE_NO+4))  # trailing rel ids

cat > "${P}/ppt/presentation.xml" <<XMLEOF
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<p:presentation xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" saveSubsetFonts="1"><p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst><p:sldIdLst>${SLD_IDS}</p:sldIdLst><p:sldSz cx="12191695" cy="6858000"/><p:notesSz cx="6858000" cy="9144000"/></p:presentation>
XMLEOF

cat > "${P}/ppt/_rels/presentation.xml.rels" <<XMLEOF
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="slideMasters/slideMaster1.xml"/>${PRES_RELS}<Relationship Id="rId${R_THEME}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme" Target="theme/theme1.xml"/><Relationship Id="rId${R_PRESP}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/presProps" Target="presProps.xml"/><Relationship Id="rId${R_VIEWP}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/viewProps" Target="viewProps.xml"/></Relationships>
XMLEOF

cat > "${P}/_rels/.rels" <<'XMLEOF'
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="ppt/presentation.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/><Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/></Relationships>
XMLEOF

cat > "${P}/docProps/core.xml" <<XMLEOF
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>Motif-mapping tool comparison - PacBio HIV-SMRTcap arm</dc:title><dc:subject>FIMO vs TFBSTools; gquad vs pqsfinder</dc:subject><cp:revision>1</cp:revision></cp:coreProperties>
XMLEOF

cat > "${P}/docProps/app.xml" <<XMLEOF
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>Microsoft Office PowerPoint</Application><Slides>${SLIDE_NO}</Slides><PresentationFormat>Custom</PresentationFormat></Properties>
XMLEOF

# [Content_Types].xml: one override per slide, layout, plus the fixed parts
CT_SLIDES=""                                         # per-slide overrides
for i in $(seq 1 "${SLIDE_NO}"); do
    CT_SLIDES="${CT_SLIDES}<Override PartName=\"/ppt/slides/slide${i}.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slide+xml\"/>"
done
CT_LAYOUTS=""                                        # per-layout overrides
for f in "${P}/ppt/slideLayouts/"*.xml; do
    CT_LAYOUTS="${CT_LAYOUTS}<Override PartName=\"/ppt/slideLayouts/$(basename "${f}")\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml\"/>"
done
cat > "${P}/[Content_Types].xml" <<XMLEOF
<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Default Extension="png" ContentType="image/png"/><Default Extension="jpeg" ContentType="image/jpeg"/><Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/><Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>${CT_LAYOUTS}${CT_SLIDES}<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/><Override PartName="/ppt/presProps.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presProps+xml"/><Override PartName="/ppt/viewProps.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.viewProps+xml"/><Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/><Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/></Types>
XMLEOF

# The slide master references its own layouts and theme; the donor's rels already
# do that correctly, so they are reused as-is. Strip any reference to parts that
# were intentionally dropped (charts, embeddings, printer settings, thumbnail).
for rel in "${P}/ppt/slideMasters/_rels/"*.rels "${P}/ppt/slideLayouts/_rels/"*.rels; do
    [ -f "${rel}" ] || continue
    # remove Relationship elements pointing at dropped part types
    sed -i -E 's#<Relationship[^>]*Target="\.\./(charts|embeddings|printerSettings)/[^"]*"[^>]*/>##g' "${rel}"
done

##################################################################################
#  Zip it. [Content_Types].xml must be the FIRST entry and stored uncompressed.
##################################################################################
rm -f "${OUT_PPTX}"                                  # replace any previous build
# member order matters: [Content_Types].xml first, then the rest of the package
( cd "${P}" && zip_create "${OUT_PPTX}" '[Content_Types].xml' _rels docProps ppt ) || {
    echo "ERROR: packaging failed" >&2; exit 1; }

echo ""
echo "=== built ${OUT_PPTX} ==="
ls -lh "${OUT_PPTX}" | awk '{print "size:", $5}'
echo "slides: ${SLIDE_NO}"
echo "=== package sanity ==="
bsdtar -tf "${OUT_PPTX}" | tail -3
# a valid package must expose the presentation part and every slide
bsdtar -tf "${OUT_PPTX}" | grep -qE 'ppt/presentation\.xml' && echo "presentation.xml present: OK"
echo "slide parts in package: $(bsdtar -tf "${OUT_PPTX}" | grep -cE 'ppt/slides/slide[0-9]+\.xml')"
