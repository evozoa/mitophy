/* site.js — tabs, hypothesis switcher, sortable tables, theme toggle */
(function () {
  'use strict';
  // theme
  const saved = localStorage.getItem('mitophy-theme');
  if (saved) document.documentElement.setAttribute('data-theme', saved);
  document.addEventListener('click', e => {
    const b = e.target.closest('[data-theme-toggle]'); if (!b) return;
    const cur = document.documentElement.getAttribute('data-theme');
    const next = cur === 'dark' ? 'light' : cur === 'light' ? '' : (matchMedia('(prefers-color-scheme: dark)').matches ? 'light' : 'dark');
    if (next) { document.documentElement.setAttribute('data-theme', next); localStorage.setItem('mitophy-theme', next); }
    else { document.documentElement.removeAttribute('data-theme'); localStorage.removeItem('mitophy-theme'); }
    document.querySelectorAll('.mtree').forEach(m => m.dispatchEvent(new Event('themechange')));
    window.dispatchEvent(new Event('resize'));
  });
  // tabs
  document.querySelectorAll('.tabs').forEach(tabs => {
    const buttons = [...tabs.querySelectorAll('button')];
    const panels = buttons.map(b => document.getElementById(b.dataset.tab));
    function activate(i, push) {
      buttons.forEach((b, j) => b.classList.toggle('active', i === j));
      panels.forEach((p, j) => p && p.classList.toggle('active', i === j));
      if (push) history.replaceState(null, '', '#' + buttons[i].dataset.tab);
      window.dispatchEvent(new Event('resize'));
      buttons[i].dispatchEvent(new CustomEvent('tabshown', { bubbles: true }));
    }
    buttons.forEach((b, i) => b.addEventListener('click', () => activate(i, true)));
    const h = location.hash.slice(1); const idx = buttons.findIndex(b => b.dataset.tab === h);
    activate(idx >= 0 ? idx : 0, false);
  });
  // sortable tables
  document.querySelectorAll('table.data.sortable').forEach(t => {
    const ths = [...t.tHead.rows[0].cells];
    ths.forEach((th, i) => th.addEventListener('click', () => {
      const asc = !(th.classList.contains('sorted') && th.classList.contains('asc'));
      ths.forEach(x => x.classList.remove('sorted', 'asc')); th.classList.add('sorted'); if (asc) th.classList.add('asc');
      const rows = [...t.tBodies[0].rows];
      const num = rows.every(r => r.cells[i] && (r.cells[i].textContent.trim() === '' || !isNaN(parseFloat(r.cells[i].textContent.replace(/,/g, '')))));
      rows.sort((a, b) => { let x = a.cells[i].textContent.trim(), y = b.cells[i].textContent.trim();
        if (num) { x = parseFloat(x.replace(/,/g, '')) || -Infinity; y = parseFloat(y.replace(/,/g, '')) || -Infinity; return asc ? x - y : y - x; }
        return asc ? x.localeCompare(y) : y.localeCompare(x); });
      rows.forEach(r => t.tBodies[0].appendChild(r));
    }));
  });
  // table filter
  document.querySelectorAll('[data-filter-table]').forEach(inp => {
    const t = document.getElementById(inp.dataset.filterTable);
    inp.addEventListener('input', () => { const q = inp.value.toLowerCase(); [...t.tBodies[0].rows].forEach(r => r.style.display = r.textContent.toLowerCase().includes(q) ? '' : 'none'); });
  });
  // hypothesis switcher: buttons[data-hyp] show .hyp-panel[data-hyp]
  document.querySelectorAll('.hyp-select').forEach(sel => {
    const btns = [...sel.querySelectorAll('button')];
    const scope = sel.closest('[data-hyp-scope]') || document;
    function show(id) { btns.forEach(b => b.classList.toggle('active', b.dataset.hyp === id));
      scope.querySelectorAll('.hyp-panel').forEach(p => { const on = p.dataset.hyp === id; p.style.display = on ? '' : 'none'; if (on && p.dataset.tree && !p.dataset.loaded) { p.dataset.loaded = 1; MTree.load(p.querySelector('.tree-slot'), p.dataset.tree, { title: p.dataset.title || '', filename: p.dataset.hyp }); } });
      window.dispatchEvent(new Event('resize')); }
    btns.forEach(b => b.addEventListener('click', () => show(b.dataset.hyp)));
    if (btns.length) show(btns[0].dataset.hyp);
  });
  // lazy trees: <div class="tree-slot" data-tree="url" data-title="..." data-collapse="1">
  document.querySelectorAll('.tree-slot[data-tree]').forEach(slot => {
    if (slot.closest('.hyp-panel')) return; // handled by switcher
    const load = () => { if (slot.dataset.loaded) return; slot.dataset.loaded = 1;
      MTree.load(slot, slot.dataset.tree, { title: slot.dataset.title || '', collapseGroups: slot.dataset.collapse === '1', filename: (slot.dataset.tree.split('/').pop() || 'tree').replace('.json', ''), supportLabels: slot.dataset.labels === '1' }); };
    const panel = slot.closest('.tabpanel');
    if (!panel || panel.classList.contains('active')) load(); else document.addEventListener('tabshown', e => { if (e.target.dataset.tab === panel.id) load(); });
  });
})();
