"""README and mkdocs pages must stay correct: their code runs, names exist, links and options are real."""
import inspect
import re
import sys
from pathlib import Path

import pytest

import nongauss
from nongauss import plots, interactive
from nongauss.cli import build_parser, main

ROOT = Path(__file__).resolve().parents[2]
PAGES = [ROOT / 'README.md', *sorted((ROOT / 'mkdocs' / 'docs').glob('*.md'))]
pytestmark = pytest.mark.skipif(not (ROOT / 'README.md').exists(), reason='needs a source checkout')


def blocks(path, lang):
    return re.findall(rf'```{lang}\n(.*?)```', path.read_text(), re.S)


@pytest.mark.parametrize('page', PAGES, ids=lambda p: p.name)
def test_python_blocks_run_in_order(page, tmp_path, monkeypatch):
    """Each page's python blocks run top to bottom in one namespace, as a reader would run them."""
    monkeypatch.chdir(tmp_path)
    ns = {}
    for code in blocks(page, 'python'):
        exec(compile(code.replace('.show()', ''), f'{page.name}', 'exec'), ns)


@pytest.mark.parametrize('page', PAGES, ids=lambda p: p.name)
def test_documented_commands_run(page, tmp_path, monkeypatch):
    """Every `nongauss ...` command shown in the docs runs (with a small grid, for speed)."""
    monkeypatch.chdir(tmp_path)
    for code in blocks(page, 'bash'):
        for line in code.splitlines():
            line = line.split('#')[0].strip()
            if not line.startswith('nongauss') or '[' in line:
                continue
            argv = line.split()[1:]
            if '--help' not in argv:
                argv += ['--N', '48']
            monkeypatch.setattr(sys, 'argv', ['nongauss', *argv])
            try:
                main()
            except SystemExit as e:             # --help exits with 0
                assert e.code == 0


def test_documented_cli_options_exist():
    known = {o for a in build_parser()._actions for o in a.option_strings}
    text = '\n'.join(p.read_text() for p in PAGES)
    cli_lines = [l for l in text.splitlines() if re.match(r'\s*(nongauss|python -m nongauss)\b', l)]
    used = set(re.findall(r'--[A-Za-z]+', ' '.join(cli_lines)))
    assert used and used <= known, used - known


def test_documented_names_exist():
    text = '\n'.join(p.read_text() for p in PAGES)
    for group in re.findall(r'from nongauss import ([\w, ]+)', text):
        for name in (n.strip() for n in group.split(',')):
            assert hasattr(nongauss, name), name
    for name in set(re.findall(r'`(plot_\w+)\(', text)):
        assert hasattr(plots, name), name
    for kind in re.findall(r'plots\.plot\(res, "(\w+)"', text):
        assert kind in plots.KINDS, kind
    assert hasattr(interactive, 'explorer')


def test_parameter_table_matches_signature():
    """Defaults listed in the Guide's generate_surface table equal the real defaults."""
    sig = inspect.signature(nongauss.generate_surface).parameters
    text = (ROOT / 'mkdocs/docs/usage.md').read_text()
    literal = {'None': None, 'True': True, 'False': False}
    checked = 0
    for names, default in re.findall(r'^\| `([\w`, ]+)` \| ([^|]+) \|', text, re.M):
        value = default.strip().strip('`')
        for n in re.findall(r'\w+', names):
            if n in sig:
                want = literal[value] if value in literal else float(value)
                assert want == sig[n].default, (n, value)
                checked += 1
    assert checked >= 9


def test_result_fields_documented():
    fields = set(nongauss.SurfaceResult.__dataclass_fields__)
    text = (ROOT / 'mkdocs/docs/usage.md').read_text()
    assert all(f'`{f}`' in text for f in fields)


def _local(target):
    """Map an image/link target of the README or docs to a path in the repo (None if external)."""
    m = re.search(r'(?:raw\.githubusercontent\.com|github\.com)/[^/]+/[^/]+/(?:blob/)?main/(\S+)', target)
    return ROOT / m.group(1) if m else None


@pytest.mark.parametrize('page', PAGES, ids=lambda p: p.name)
def test_images_and_links_exist(page):
    text = page.read_text()
    for target in re.findall(r'!\[[^\]]*\]\(([^)\s]+)', text) + re.findall(r'src="([^"]+)"', text) \
            + re.findall(r'\]\((https://htmlpreview[^)]+)\)', text):
        local = _local(target)
        if local is None and not target.startswith('http'):
            local = (page.parent / target).resolve()      # relative path inside the docs
        if local is not None:
            assert local.exists(), f'{page.name}: {target}'


def test_nav_pages_exist():
    cfg = (ROOT / 'mkdocs' / 'config.yml').read_text()
    pages = re.findall(r':\s+(\S+\.md)\s*$', cfg, re.M)
    assert pages and all((ROOT / 'mkdocs' / 'docs' / p).exists() for p in pages)
