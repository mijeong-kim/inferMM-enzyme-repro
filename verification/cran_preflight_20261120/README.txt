CRAN preflight and PDF rendering record — 20 November 2026

The current tar is identified in packages/package_build_current.csv. The
final_report.json and source_diff.json are verbatim records from the authors'
local preflight. Their absolute paths identify that original execution.
This package has 0 check errors, 0 warnings and 1 local environment NOTE:
HTML validation was skipped because recent HTML Tidy and R package V8 were
unavailable. Online CRAN incoming, examples, vignettes, PDF and 694 tests passed.
No CRAN submission or remote check upload was performed.

The source keeps canonical Rd DOI markup. On the authors' local R/TeX,
discretionary hyphens in DOI display labels disappeared from the stock PDF.
The helper restores only visible hyphen tokens in transient LaTeX labels;
hyperlink destinations and package sources stay unchanged. The distributed
helper adds --source/--output/--qa arguments to the reviewed local renderer.
It requires Python 3, R and TeX (pdflatex and makeindex), not network access.

To reproduce the supplied manual from this extracted bundle, use a scratch
directory and unpack packages/inferMM_0.1.2.tar.gz there. For example:
  mkdir /tmp/infermm-manual-source
  tar -xzf packages/inferMM_0.1.2.tar.gz -C /tmp/infermm-manual-source
  python3 verification/cran_preflight_20261120/render_manual_preserve_doi_hyphens.py --source /tmp/infermm-manual-source/inferMM --output /tmp/inferMM_0.1.2_manual.pdf --qa /tmp/infermm-manual-qa
PDF creation timestamps and platform fonts can change PDF bytes; the checks
record the supplied PDF's 26 pages, DOI labels, links and visual review.
Large page/contact PNGs and local build intermediates are not in this ZIP.

The numerical estimation R files, tests, NAMESPACE, demonstration data, analysis
sources, original 0.1.1 execution records, saved results and timing are retained.
Earlier build-specific validation fields and source hashes remain historical.
