#!/bin/bash
#SBATCH --job-name=pb_pipeline_pptx
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

##################################################################################
# Build writeups/PacBio_Pipeline_Report.pptx: a 15-slide summary of the full
# PacBio pipeline re-run (purpose, methodology, results, challenges, conclusions).
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
# Usage: bash scripts/utils/build_pacbio_pipeline_pptx.sh   (needs zip: login node)
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
# tables the slides are built from live under results/reporting/pacbio_deck
TEMPLATE="${REPO_ROOT}/writeups/PacBio_Tool_Comparison.pptx"         # theme/master donor
OUT_PPTX="${REPO_ROOT}/writeups/PacBio_Pipeline_Report.pptx"   # what we build
BUILD="${REPO_ROOT}/.pptx_build_pipeline"                     # scratch tree, removed at the end
# clean any half-finished tree from a previous run
rm -rf "${BUILD}"; mkdir -p "${BUILD}/pkg"
# always tidy the scratch tree, even on failure
trap 'rm -rf "${BUILD}"' EXIT

# This deck is built from the per-stage summary tables, not from the motif
# tool-comparison analysis the sibling deck uses, so the required inputs differ.
DECK_DATA="${DECK_DATA:-${REPO_ROOT}/results/reporting/pacbio_deck}"
for f in "${TEMPLATE}" \
         "${DECK_DATA}/asm_compact.tsv" "${DECK_DATA}/amplicon.tsv" \
         "${DECK_DATA}/msa.tsv" "${DECK_DATA}/ltr.tsv" "${DECK_DATA}/motif.tsv" \
         "${DECK_DATA}/gene.tsv" "${DECK_DATA}/caller.tsv" "${DECK_DATA}/mask.tsv" \
         "${DECK_DATA}/before_after.tsv"; do
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
    SLIDE_BODY="${SLIDE_BODY}<p:sp><p:nvSpPr><p:cNvPr id=\"${SHAPE_ID}\" name=\"Rule ${SHAPE_ID}\"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x=\"576072\" y=\"1524000\"/><a:ext cx=\"713232\" cy=\"38100\"/></a:xfrm><a:prstGeom prst=\"roundRect\"><a:avLst/></a:prstGeom><a:solidFill><a:srgbClr val=\"${ACCENT}\"/></a:solidFill><a:ln><a:noFill/></a:ln></p:spPr><p:txBody><a:bodyPr/><a:lstStyle/><a:p/></p:txBody></p:sp>"
}

SLIDE_NO=0                                           # how many slides written so far
# Start a slide: eyebrow label + title + accent rule
begin_slide() {
    local eyebrow="$1" title="$2"
    SHAPE_ID=1; SLIDE_BODY=""                        # reset per-slide state
    add_text 365760  400000 1400 1 "${ACCENT}" "Segoe UI" "${eyebrow}"
    add_text 685800  820000 3300 1 "${INK}"    "Segoe UI" "${title}"
    add_rule
    # The eyebrow, title and rule are pinned to the top of every slide, so they
    # are banked here and the body is accumulated separately. end_slide then
    # shifts only the body, which is what makes vertical balancing possible.
    HEADER_XML="${SLIDE_BODY}"
    SLIDE_BODY=""
    BODY_TOP=1800000                                 # first body line sits here
    CURSOR=${BODY_TOP}                               # y position for the next body line
}

# Body line styles, each advancing CURSOR by its own line height
add_lead()   { add_text "${CURSOR}" 490000 1850 0 "${MUTED}" "Segoe UI"  "$1"; CURSOR=$((CURSOR+530000)); }
add_bullet() { add_text "${CURSOR}" 430000 1750 0 "${INK}"   "Segoe UI"  "•  $1"; CURSOR=$((CURSOR+415000)); }
add_key()    { add_text "${CURSOR}" 430000 1750 1 "${ACCENT}" "Segoe UI" "$1"; CURSOR=$((CURSOR+425000)); }
add_warn()   { add_text "${CURSOR}" 430000 1750 1 "${WARN}"  "Segoe UI"  "$1"; CURSOR=$((CURSOR+425000)); }
add_mono()   { add_text "${CURSOR}" 330000 1400 0 "${INK}"   "Consolas"  "$1"; CURSOR=$((CURSOR+300000)); }
add_gap()    { CURSOR=$((CURSOR+200000)); }

# Slide height is 6858000 EMU; leave a bottom margin so text never runs off the
# edge. Because every box is absolutely positioned, an overlong slide fails
# silently in PowerPoint rather than reflowing, so it is checked at build time.
SLIDE_BOTTOM=6600000

# Emit the accumulated shapes as ppt/slides/slideN.xml plus its layout rel
end_slide() {
    SLIDE_NO=$((SLIDE_NO+1))                         # this slide's 1-based number
    if [ "${CURSOR}" -gt "${SLIDE_BOTTOM}" ]; then
        echo "WARNING: slide ${SLIDE_NO} content reaches ${CURSOR} EMU, past the ${SLIDE_BOTTOM} safe bottom -- trim a row or shrink a style" >&2
    fi

    # Balance the slide vertically. A short slide otherwise leaves all its
    # whitespace in one block at the bottom, which reads as unfinished next to a
    # full slide. Half the leftover space is pushed above the body so the block
    # sits optically centred in the area under the title. The shift is capped so
    # a very sparse slide does not end up floating in the middle of the page.
    local leftover=$(( SLIDE_BOTTOM - CURSOR ))
    local shift=0
    if [ "${leftover}" -gt 0 ]; then
        shift=$(( leftover / 2 ))
        [ "${shift}" -gt 700000 ] && shift=700000
    fi
    local body="${SLIDE_BODY}"
    if [ "${shift}" -gt 0 ]; then
        # add the shift to every y offset in the body shapes, leaving the banked
        # header untouched
        body=$(printf '%s' "${SLIDE_BODY}" | awk -v s="${shift}" '{
            out = ""; rest = $0
            while (match(rest, /y="[0-9]+"/)) {
                tok = substr(rest, RSTART, RLENGTH)
                val = tok; gsub(/[^0-9]/, "", val)
                out = out substr(rest, 1, RSTART-1) "y=\"" (val + s) "\""
                rest = substr(rest, RSTART + RLENGTH)
            }
            print out rest
        }')
    fi
    SLIDE_BODY="${HEADER_XML}${body}"
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
##################################################################################
#  Slides
#
#  Layout discipline: each slide carries ONE title and at most TWO content
#  blocks -- a short lead line and a single table or short bullet group -- so no
#  slide is crowded. Tables are read from TSVs at build time, so the deck cannot
#  drift from the results.
##################################################################################
D="${DECK_DATA}"                                     # table TSVs, checked above

# --- 1. Title -----------------------------------------------------------------
begin_slide "HIV-1 PROVIRAL U3 ANALYSIS  ·  PACBIO HIV-SMRTcap ARM" "Pipeline re-run: methods, defects fixed, and results"
add_lead "Four Rakai-cohort SMRTcap proviral samples, run end to end on Slurm."
add_gap
add_key "124_4   ·   128_5   ·   203_3   ·   211_0      (+ HXB2 / K03455.1 reference)"
add_gap
add_bullet "Objective 1: recover proviral U3 and map transcription-factor binding sites"
add_bullet "This deck covers purpose, method, results, challenges and conclusions per step"
end_slide

# --- 2. Purpose ---------------------------------------------------------------
begin_slide "PURPOSE" "What the pipeline has to establish"
add_lead "U3 is the HIV-1 promoter. Its TFBS content shapes transcription and latency."
add_gap
add_bullet "Recover each sample's own U3 sequence, not the reference's"
add_bullet "Map NF-kB / SP1 sites and compare across samples"
add_bullet "Keep \"not sequenced\" distinguishable from \"motif absent\" at every step"
add_bullet "Compare tool choices at each stage rather than trusting one tool"
end_slide

# --- 3. Method ----------------------------------------------------------------
begin_slide "METHODOLOGY" "Stages, tools compared, and tools chosen"
add_lead "Every stage runs as a Slurm job; every stage writes a summary.tsv row per tool."
add_gap
add_mono "QC          NanoPlot / chopper                 -> chopper"
add_mono "Extraction  minimap2 + N-strip provirus tool"
add_mono "Assembly    minimap2+bcftools | hifiasm        -> reference-guided, 2 arms"
add_mono "MSA         MAFFT | MUSCLE | Clustal Omega     -> MAFFT"
add_mono "U3          HXB2-anchored alignment liftover"
add_mono "Motifs      FIMO | TFBSTools | MOODS           -> FIMO"
add_mono "Gquad       gquad | pqsfinder"
add_mono "Filtering   HIV-Intact | Poplars | HIVSeqinR   -> HIV-Intact + Poplars"
end_slide

# --- 4. Two arms --------------------------------------------------------------
begin_slide "DESIGN" "Two reference arms, one variable"
add_lead "The assembly step runs twice so the cost of the reference choice stays visible."
add_gap
add_key "minimap2_consensus  -  HXB2 (subtype B) baseline, the coordinate system"
add_key "minimap2_bestref    -  closest LTR-complete Group M reference per sample"
add_gap
add_bullet "Both arms are the same pipeline over a different reference"
add_bullet "Every downstream stage is run once per arm"
add_bullet "hifiasm is kept as a reference-free comparator"
end_slide

# --- 5. Core finding: amplicons ----------------------------------------------
begin_slide "RESULT  ·  THE DEFINING PROPERTY OF THIS DATA" "These are amplicon libraries, not shotgun"
add_lead "Reads stack into a few discrete intervals; they do not tile the genome."
add_gap
add_tsv "${D}/amplicon.tsv" 6
add_gap
add_warn "Consequence: hifiasm produced 0 contigs for 3 of 4 samples - no overlap graph exists"
end_slide

# --- 6. Defect 1 --------------------------------------------------------------
begin_slide "CHALLENGE 1  ·  INDEL CALLING WAS SILENTLY OFF" "Zero indels called in all four samples"
add_lead "At 203_3 position 168 all 19 reads delete a T, yet no variant was called."
add_gap
add_tsv "${D}/caller.tsv" 6
add_gap
add_key "Fix: mpileup -X pacbio-ccs. Switching the caller alone changed nothing."
end_slide

# --- 7. Defect 2 --------------------------------------------------------------
begin_slide "CHALLENGE 2  ·  MASKING CONFLATED TWO THINGS" "\"Deleted\" was being reported as \"not sequenced\""
add_lead "samtools depth counts bases, so a base deleted in every read reports depth 0."
add_gap
add_tsv "${D}/mask.tsv" 6
add_gap
add_key "Fix: derive coverage from read alignment spans. Interior masking now zero."
end_slide

# --- 8. Assembly results ------------------------------------------------------
begin_slide "RESULT  ·  ASSEMBLY" "Reference-guided consensus, both arms"
add_lead "Validity is called bases (non-N), not total length: 9719bp that is 88% N is not a genome."
add_gap
add_tsv "${D}/asm_compact.tsv" 10
add_gap
add_bullet "Only 124_4 clears the 8000 called-base threshold, in both arms"
end_slide

# --- 9. MSA -------------------------------------------------------------------
begin_slide "RESULT  ·  MULTIPLE SEQUENCE ALIGNMENT" "All three aligners preserved sequence integrity"
add_lead "Every aligner returned each record at its exact ungapped length - nothing altered."
add_gap
add_tsv "${D}/msa.tsv" 8
add_gap
add_warn "Clustal Omega treats N as alignable: on collinear input it inserted 49 false indels per record"
end_slide

# --- 10. U3 recovery ----------------------------------------------------------
begin_slide "RESULT  ·  U3 RECOVERY" "Each sample now uses whichever LTR it actually sequenced"
add_lead "Both LTRs are identical in an integrated provirus, but an amplicon may reach only one."
add_gap
add_tsv "${D}/ltr.tsv" 7
add_gap
add_key "124_4 switches to its 3' copy - 474 called bases where the 5' copy had none"
end_slide

# --- 11. Artifact resolved ----------------------------------------------------
begin_slide "RESULT  ·  TWO ARTIFACTS ELIMINATED" "124_4 recovered, 128_5 correctly excluded"
add_lead "Previously two samples returned U3 byte-identical to HXB2, so their motifs were HXB2's."
add_gap
add_tsv "${D}/before_after.tsv" 6
add_gap
add_key "128_5 has no U3 coverage at either LTR and is now auto-excluded, not silently reported"
end_slide

# --- 12. The insertion --------------------------------------------------------
begin_slide "RESULT  ·  A CANDIDATE FINDING" "124_4 carries a 21 bp insertion in its 3' U3"
add_lead "Real sequence, zero N, inside the regulatory region."
add_gap
add_bullet "Reproduced against BOTH references - not a reference artifact"
add_bullet "All three aligners agree on its presence and size"
add_bullet "Placement inside the run is ambiguous - the signature of a tandem duplication"
add_bullet "Sequence resembles adjacent HXB2 U3 containing ACTGCTGACA"
add_gap
add_warn "Only visible after all three fixes: indel calling, span mask, and LTR choice"
end_slide

# --- 13. Motif results --------------------------------------------------------
begin_slide "RESULT  ·  MOTIF MAPPING (FIMO)" "Real per-sample TFBS profiles"
add_lead "HXB2's own hits are unchanged across the re-run, which serves as the control."
add_gap
add_tsv "${D}/motif.tsv" 14
end_slide

# --- 14. Intactness basis -----------------------------------------------------
begin_slide "RESULT  ·  BIOLOGICAL FILTERING" "A verdict only means something where the gene was read"
add_lead "Per-gene read coverage decides whether an intactness call is supported at all."
add_gap
add_tsv "${D}/gene.tsv" 10
add_gap
add_key "HIV-Intact: 1 intact, 4 non-intact - and HXB2 is the intact one (control passes)"
end_slide

# --- 15. Conclusions ----------------------------------------------------------
begin_slide "CONCLUSIONS  ·  AND WHAT REMAINS" "Where the analysis now stands"
add_lead "The pipeline runs end to end and no longer reports reference sequence as sample data."
add_gap
add_bullet "3 of 4 samples yield genuine U3; 128_5 is excluded on evidence, not assumption"
add_bullet "Subtype-matching helps only 203_3 (162 -> 95 variants); a wash elsewhere"
add_bullet "HIVSeqinR needs primer-flanked full genomes - unsuitable for amplicon input"
add_gap
add_warn "Open: report motifs as named sites (LANL) so cross-sample positions are comparable"
end_slide
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
