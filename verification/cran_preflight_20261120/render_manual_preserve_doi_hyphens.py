"""Render the canonical Rd manual, preserving DOI display hyphens in local TeX.

The local R tools:::Rd_expr_doi expansion uses TeX discretionary hyphens
(\\-) in DOI labels. They disappear when no line break occurs. This helper
changes only those display-label tokens in transient LaTeX. Canonical package
Rd files, DOI hyperlink targets, R namespaces, and the R installation stay intact.
Run with Python 3 from any working directory; R and a TeX installation are needed.
"""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import tempfile

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, required=True, help='Unpacked inferMM source directory')
parser.add_argument('--output', type=Path, default=Path('inferMM_0.1.2_manual.pdf'))
parser.add_argument('--qa', type=Path, default=None)
args = parser.parse_args()
output_target = args.output.resolve()
PREFLIGHT = output_target.parent
PREFLIGHT.mkdir(parents=True, exist_ok=True)
SOURCE = args.source.resolve()
assert (SOURCE / 'DESCRIPTION').is_file(), SOURCE
QA = args.qa.resolve() if args.qa else PREFLIGHT / 'pdf_qa'
QA.mkdir(parents=True, exist_ok=True)
QA.mkdir(exist_ok=True)
BUILD = Path(tempfile.mkdtemp(prefix='doi_latex_', dir=QA))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
source_hashes = {str(p.relative_to(SOURCE)): sha(p)
                 for p in SOURCE.rglob('*') if p.is_file()}
stock_pdf = QA / 'manual_stock_from_canonical_Rd.pdf'
with (QA / 'doi_stock_build.log').open('w') as log:
    subprocess.run(['R', 'CMD', 'Rd2pdf', '--no-clean', '--no-preview', '--force',
                    '--build-dir=' + str(BUILD), '--output=' + str(stock_pdf),
                    str(SOURCE)], cwd=PREFLIGHT, stdout=log,
                   stderr=subprocess.STDOUT, check=True)
tex = BUILD / 'Rd2.tex'
before = tex.read_text()
(BUILD / 'Rd2.canonical.tex').write_text(before)

def braced(text, start):
    assert text[start] == '{'
    depth = 1
    pos = start + 1
    while depth:
        depth += (text[pos] == '{') - (text[pos] == '}')
        pos += 1
    return text[start + 1:pos - 1], pos

replacements = []
targets = []
pos = 0
while True:
    pos = before.find('\\Rhref{', pos)
    if pos < 0:
        break
    target, next_group = braced(before, pos + len('\\Rhref'))
    label, end = braced(before, next_group)
    if target.startswith('https://doi.org/') and label.startswith('doi:'):
        targets.append(target)
        patched = label.replace('\\-', '-')
        if patched != label:
            replacements.append((next_group + 1, end - 1, label, patched, target))
    pos = end
assert len(targets) == 13, len(targets)
assert len(replacements) == 6, len(replacements)
assert sum(r[2].count('\\-') for r in replacements) == 14
after = before
for start, end, old, new, target in reversed(replacements):
    assert after[start:end] == old
    after = after[:start] + new + after[end:]
tex.write_text(after)

# Match R's per-process Rd style lookup without touching an installed style file.
rhome = subprocess.check_output(['R', 'RHOME'], text=True).strip()
env = os.environ.copy()
env['TEXINPUTS'] = str(Path(rhome) / 'share/texmf/tex/latex') + '//' + os.pathsep + env.get('TEXINPUTS', '')
pdflatex = shutil.which('pdflatex')
makeindex = shutil.which('makeindex')
assert pdflatex and makeindex
commands = [[pdflatex, '-interaction=nonstopmode', '-halt-on-error', 'Rd2.tex'],
            [makeindex, 'Rd2.idx'],
            [pdflatex, '-interaction=nonstopmode', '-halt-on-error', 'Rd2.tex'],
            [pdflatex, '-interaction=nonstopmode', '-halt-on-error', 'Rd2.tex']]
with (QA / 'doi_corrected_latex_build.log').open('w') as log:
    for command in commands:
        subprocess.run(command, cwd=BUILD, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
assert all(sha(SOURCE / rel) == expected for rel, expected in source_hashes.items())
output = output_target
staged = QA / 'manual_corrected_staged.pdf'
shutil.copy2(BUILD / 'Rd2.pdf', staged)
staged.replace(output)
report = dict(status='RENDERED_REQUIRES_VISUAL_QA', output=str(output),
              pdf_sha256=sha(output), build_directory=str(BUILD),
              canonical_Rd_unchanged=True, all_package_source_files_unchanged=True,
              global_R_namespace_or_installation_changed=False,
              DOI_hyperlink_targets_preserved=targets,
              corrected_label_count=len(replacements), corrected_hyphens=14,
              labels=[dict(target=r[4], before=r[2], after=r[3]) for r in replacements])
(QA / 'doi_render_patch.json').write_text(json.dumps(report, ensure_ascii=False, indent=2))
print(json.dumps({k:report[k] for k in ['status', 'output', 'pdf_sha256', 'corrected_label_count',
                                      'corrected_hyphens', 'all_package_source_files_unchanged']}, indent=2))
