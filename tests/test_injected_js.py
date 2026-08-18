#!/usr/bin/env python
"""Execute the app's injected browser scripts the way the browser does, and check they survive.

    ~/.conda/envs/cluster_annotation/bin/python tests/test_injected_js.py

`st.html` injects inline, so every script the app injects shares ONE global lexical scope — and is
re-injected on every rerun. Two scripts that both declare `const doc` at top level, or one script
run twice, is a *parse-time* SyntaxError: the script dies whole, listeners and all, with nothing in
the UI to say so. That is precisely how the keyboard shortcuts broke — they registered zero
listeners while every other part of the app looked healthy.

So the scripts are pulled out of app.py as the app emits them, wrapped exactly as `inject_js`
wraps them, and run in one shared context under node. Requires `node`; skips without it.
"""
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
APP = ROOT / "app" / "app.py"
HAVE_NODE = shutil.which("node") is not None

# A minimal DOM: enough for the scripts to register listeners and find their targets.
HARNESS = r"""
const vm = require('vm'), fs = require('fs');
const scripts = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const listeners = [];
const badge = { textContent: 'not detected' };
const clicked = [];
const doc = {
  addEventListener: (t, f, c) => listeners.push({type: t, fn: f, capture: !!c}),
  removeEventListener: () => {},
  getElementById: (id) => (id === 'kb-status' ? badge : null),
  querySelector: (sel) => ({ disabled: false, click: () => clicked.push(sel) }),
};
const win = { document: doc, getSelection: () => '', requestAnimationFrame: (f) => f() };
win.parent = win;
const ctx = vm.createContext({ window: win, document: doc, console });
const errors = [];
// Run every script twice: a rerun re-injects each one into the same scope.
for (const pass of [1, 2]) {
  scripts.forEach((s, i) => {
    try { vm.runInContext(s, ctx); }
    catch (e) { errors.push(`pass ${pass} script ${i + 1}: ${e.name}: ${e.message}`); }
  });
}
// Fire the shortcuts through whichever keydown handlers registered.
const press = (key) => {
  const ev = { key, target: { tagName: 'BODY' }, preventDefault: () => {}, stopImmediatePropagation: () => {} };
  listeners.filter(l => l.type === 'keydown').forEach(l => l.fn(ev));
};
['ArrowLeft', 'ArrowRight', '[', ']', 'u', '1', '0'].forEach(press);
console.log(JSON.stringify({
  errors, listeners: listeners.length, clicked, badge: badge.textContent,
}));
"""


def _scripts() -> list[str]:
    """The scripts app.py injects, rendered and wrapped exactly as inject_js does."""
    src = APP.read_text()
    # Mirror whatever isolation inject_js applies, rather than asserting one exact spelling: if it
    # stops isolating, the scripts below fail behaviourally (a real SyntaxError), which is the
    # failure worth reporting.
    inject = re.search(r"def inject_js.*?\n\n\ndef ", src, re.S).group(0)
    isolates = "(function(){" in inject or "=>{" in inject or "=> {" in inject
    bodies = re.findall(r'inject_js\((f?)"""(.*?)"""\)', src, re.S)
    assert bodies, "no injected scripts found in app.py"
    out = []
    for is_f, body in bodies:
        if is_f:      # the keymap is interpolated; any section list exercises the same code
            sections = [f"s{i}" for i in range(10)]
            keymap = {"ArrowLeft": ".st-key-fp_prev button", "ArrowRight": ".st-key-fp_next button",
                      "[": ".st-key-cl_prev button", "]": ".st-key-cl_next button",
                      "u": ".st-key-next_unrev button"}
            for i, _ in enumerate(sections):
                keymap[str((i + 1) % 10)] = f".st-key-sec_{i} button"
            body = eval('f"""' + body + '"""', {"json": json, "keymap": keymap})  # noqa: S307
        out.append("(function(){\n" + body + "\n})();" if isolates else body)
    return out


def _run() -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        sp, hp = Path(tmp) / "scripts.json", Path(tmp) / "harness.js"
        sp.write_text(json.dumps(_scripts()))
        hp.write_text(HARNESS)
        res = subprocess.run(["node", str(hp), str(sp)], capture_output=True, text=True)
        assert res.returncode == 0, res.stderr
        return json.loads(res.stdout)


def test_every_injected_script_parses_and_survives_a_rerun():
    """The regression: the second script used to die on `const doc` already being declared."""
    if not HAVE_NODE:
        print("      (no node; skipping)")
        return
    got = _run()
    assert not got["errors"], "injected script(s) failed to run:\n  " + "\n  ".join(got["errors"])


def test_the_keydown_listeners_actually_register():
    if not HAVE_NODE:
        print("      (no node; skipping)")
        return
    got = _run()
    # One listener per script per pass: the cache-key blocker and the shortcut handler.
    assert got["listeners"] >= 2, got


def test_each_shortcut_clicks_its_button():
    if not HAVE_NODE:
        print("      (no node; skipping)")
        return
    clicked = _run()["clicked"]
    for sel in (".st-key-fp_prev button", ".st-key-fp_next button", ".st-key-cl_prev button",
                ".st-key-cl_next button", ".st-key-next_unrev button", ".st-key-sec_0 button",
                ".st-key-sec_9 button"):
        assert sel in clicked, f"{sel} was never clicked; got {sorted(set(clicked))}"


def test_the_status_badge_reports_the_listener_is_live():
    """The badge is what makes a future silent failure visible in the UI."""
    if not HAVE_NODE:
        print("      (no node; skipping)")
        return
    assert "active" in _run()["badge"]


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"  ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
