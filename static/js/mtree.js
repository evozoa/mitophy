/* mtree.js — minimal dependency-free phylogenetic tree viewer (SVG).
 * Input: tree JSON {newick, tips:{label:{name,group,color,accession,...}}, nodes:{id:{ufboot,alrt,support,note,label}}, legend:[{id,color}]}
 * Usage: MTree.render(containerElement, treeJson, {phylogram:true, tipHeight:15})
 * (c) mitophy, MIT licence */
(function (global) {
  'use strict';

  // ---------- Newick parsing ----------
  function parseNewick(s) {
    let i = 0;
    const root = { children: [], length: null, label: '' };
    function readLabel() {
      let out = '';
      if (s[i] === "'") { i++; while (i < s.length) { if (s[i] === "'" && s[i + 1] === "'") { out += "'"; i += 2; } else if (s[i] === "'") { i++; break; } else out += s[i++]; } return out; }
      while (i < s.length && !'(),:;[]'.includes(s[i])) out += s[i++];
      return out.trim();
    }
    function readNode(node) {
      skipWs();
      if (s[i] === '(') {
        i++;
        while (true) {
          const child = { children: [], length: null, label: '' };
          readNode(child); node.children.push(child); skipWs();
          if (s[i] === ',') { i++; continue; }
          if (s[i] === ')') { i++; break; }
          throw new Error('Newick parse error at ' + i);
        }
      }
      skipWs(); node.label = readLabel(); skipWs();
      if (s[i] === '[') { while (i < s.length && s[i] !== ']') i++; i++; }
      if (s[i] === ':') { i++; skipWs(); let n = ''; while (i < s.length && !'(),;[]'.includes(s[i])) n += s[i++]; node.length = parseFloat(n); }
      skipWs();
      if (s[i] === '[') { while (i < s.length && s[i] !== ']') i++; i++; }
    }
    function skipWs() { while (i < s.length && /\s/.test(s[i])) i++; }
    readNode(root);
    return root;
  }

  // ---------- helpers ----------
  const SVGNS = 'http://www.w3.org/2000/svg';
  function el(tag, attrs, parent) {
    const e = document.createElementNS(SVGNS, tag);
    for (const k in attrs) if (attrs[k] !== undefined && attrs[k] !== null) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function html(tag, attrs, parent, text) {
    const e = document.createElement(tag);
    for (const k in attrs) if (k === 'class') e.className = attrs[k]; else e.setAttribute(k, attrs[k]);
    if (text !== undefined) e.textContent = text;
    if (parent) parent.appendChild(e);
    return e;
  }
  function esc(t) { return String(t).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c])); }

  function supportText(n) {
    if (!n) return '';
    if (n.ufboot !== undefined && n.alrt !== undefined) return 'SH-aLRT ' + n.alrt + ' / UFBoot ' + n.ufboot;
    if (n.support !== undefined) return 'support ' + n.support;
    if (n.label) return n.label;
    return '';
  }
  function supportValue(n) { // 0..100 or null
    if (!n) return null;
    if (n.ufboot !== undefined) return n.ufboot;
    if (n.support !== undefined) return n.support;
    return null;
  }

  // ---------- main ----------
  function render(container, data, opts) {
    opts = Object.assign({ phylogram: true, tipHeight: 15, fontSize: 11, labelWidth: 230, minWidth: 700, showSupport: true,
      supportThreshold: 0, collapseGroups: false, title: '' }, opts || {});
    container.innerHTML = '';
    container.classList.add('mtree');
    const tips = data.tips || {}, nodesMeta = data.nodes || {};
    const root = parseNewick(data.newick);
    let uid = 0;
    (function walk(n, parent, depth) { n.id = uid++; n.parent = parent; n.depth = depth; n.collapsed = false;
      n.meta = n.children.length ? (nodesMeta[n.label] || {}) : (tips[n.label] || { name: n.label });
      n.children.forEach(c => walk(c, n, depth + 1)); })(root, null, 0);
    const hasLengths = (function anyLen(n) { return (n.length !== null && !isNaN(n.length) && n.length > 0) || n.children.some(anyLen); })(root);
    let phylogram = opts.phylogram && hasLengths;
    let query = '';

    // toolbar
    const bar = html('div', { class: 'mtree-bar' }, container);
    if (opts.title) html('span', { class: 'mtree-title' }, bar, opts.title);
    const search = html('input', { type: 'search', placeholder: 'find taxon…', class: 'mtree-search' }, bar);
    const btnLen = html('button', { type: 'button', class: 'mtree-btn' }, bar, phylogram ? 'cladogram' : 'branch lengths');
    if (!hasLengths) btnLen.disabled = true;
    const btnExpand = html('button', { type: 'button', class: 'mtree-btn' }, bar, 'expand all');
    const btnCollapse = html('button', { type: 'button', class: 'mtree-btn' }, bar, 'collapse by group');
    const btnSvg = html('button', { type: 'button', class: 'mtree-btn' }, bar, 'download SVG');
    const info = html('span', { class: 'mtree-info' }, bar, '');

    // legend
    if (data.legend && data.legend.length) {
      const lg = html('div', { class: 'mtree-legend' }, container);
      data.legend.forEach(g => {
        const it = html('span', { class: 'mtree-legend-item' }, lg);
        const sw = html('span', { class: 'mtree-swatch' }, it); sw.style.background = g.color;
        it.appendChild(document.createTextNode(g.id));
        it.title = 'click to collapse/expand this group';
        it.addEventListener('click', () => { toggleGroup(g.id); draw(); });
      });
    }
    const wrap = html('div', { class: 'mtree-svgwrap' }, container);
    const tooltip = html('div', { class: 'mtree-tooltip' }, container);
    tooltip.style.display = 'none';

    function toggleGroup(gid) {
      // collapse maximal clades whose tips all belong to gid; if any such clade is collapsed, expand instead
      let any = false;
      (function scan(n) { if (!n.children.length) return; if (n.collapsed && groupOf(n) === gid) any = true; n.children.forEach(scan); })(root);
      (function walk(n) {
        if (!n.children.length) return;
        if (groupOf(n) === gid && n.parent && groupOf(n.parent) !== gid) { n.collapsed = !any; return; }
        n.children.forEach(walk);
      })(root);
    }
    function groupOf(n) { // group id if all tips share it, else null (cached)
      if (n._grp !== undefined) return n._grp;
      if (!n.children.length) return (n._grp = (n.meta.group || null));
      const gs = new Set(n.children.map(groupOf));
      return (n._grp = (gs.size === 1 ? [...gs][0] : null));
    }
    function colorOf(n) { const g = groupOf(n); if (!g) return null; if (!n.children.length) return n.meta.color || null;
      let t = n; while (t.children.length) t = t.children[0]; return t.meta.color || null; }
    function tipCount(n) { return n.children.length ? n.children.reduce((a, c) => a + tipCount(c), 0) : 1; }
    function displayName(n) { return n.meta.name || n.label; }

    function layout() {
      // visible leaves (collapsed clades count as one row)
      let row = 0; const vis = [];
      (function walk(n) {
        n.x = 0;
        if (!n.children.length || n.collapsed) { n.y = row++; vis.push(n); return; }
        n.children.forEach(walk);
        n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2;
      })(root);
      // x: distance from root
      (function walk(n, d, depthMax) {
        n.dist = d; n.children.forEach(c => walk(c, d + (phylogram ? (c.length || 0) : 1)));
      })(root, 0);
      let maxD = 0; (function walk(n) { if (!n.children.length || n.collapsed) maxD = Math.max(maxD, n.dist); n.children.forEach(walk); })(root);
      if (!phylogram) { // cladogram: align tips right
        (function walk(n) { if (!n.children.length || n.collapsed) n.dist = maxD; else n.children.forEach(walk); })(root);
        // parents at min(child)-1 for a tidy look
        (function walk(n) { if (n.children.length && !n.collapsed) { n.children.forEach(walk); n.dist = Math.min(...n.children.map(c => c.dist)) - 1; } })(root);
        const minD = root.dist; (function walk(n) { n.dist -= minD; n.children.forEach(walk); })(root); maxD -= minD;
      }
      return { rows: row, maxD: maxD || 1, vis };
    }

    function draw() {
      const L = layout();
      const th = opts.tipHeight, fs = opts.fontSize;
      const availW = Math.max(opts.minWidth, wrap.clientWidth || 800);
      const padL = 12, padT = 24, padB = 30, labelW = opts.labelWidth;
      const treeW = Math.max(200, availW - padL - labelW - 30);
      const H = L.rows * th + padT + padB, W = availW;
      const sx = d => padL + d / L.maxD * treeW, sy = r => padT + r * th + th / 2;
      wrap.innerHTML = '';
      const svg = el('svg', { xmlns: SVGNS, width: W, height: H, viewBox: `0 0 ${W} ${H}`, class: 'mtree-svg', 'font-family': 'system-ui, sans-serif', 'font-size': fs }, wrap);
      const style = el('style', {}, svg);
      style.textContent = '.mt-edge{fill:none;stroke:currentColor;stroke-width:1.2}.mt-tip{cursor:default}.mt-tip.hit text{font-weight:700;fill:#c0392b}.mt-node{cursor:pointer}.mt-node:hover circle{stroke-width:2}.mt-sup{font-size:' + (fs - 2) + 'px;fill:currentColor;opacity:.75}.mt-clade{cursor:pointer}';
      const gEdges = el('g', { class: 'mt-edges', color: 'currentColor' }, svg);
      const gNodes = el('g', {}, svg);
      const gTips = el('g', {}, svg);
      const q = query.toLowerCase();
      let hits = 0;
      (function walk(n) {
        const x = sx(n.dist), y = sy(n.y);
        if (n.parent) {
          const px = sx(n.parent.dist), py = sy(n.parent.y);
          el('path', { class: 'mt-edge', d: `M${px},${py}V${y}H${x}` }, gEdges);
        }
        if (n.children.length && !n.collapsed) {
          n.children.forEach(walk);
          // internal node: support marker + tooltip
          const sv = supportValue(n.meta), st = supportText(n.meta), note = n.meta.note;
          const g = el('g', { class: 'mt-node' }, gNodes);
          if (opts.showSupport && (sv !== null || note || n.meta.label)) {
            if (sv !== null && sv >= opts.supportThreshold) {
              const fill = sv >= 95 ? '#222' : sv >= 80 ? '#888' : '#fff';
              el('circle', { cx: x, cy: y, r: 3.2, fill: fill, stroke: '#222', 'stroke-width': 1 }, g);
            }
            // print the value only where it is informative (below 95); a filled dot already means >= 95
            if (sv !== null && sv < 95 && (opts.supportLabels || L.rows <= 60)) el('text', { x: x - 5, y: y - 3, 'text-anchor': 'end', class: 'mt-sup' }, g).textContent = (n.meta.ufboot !== undefined && n.meta.alrt !== undefined) ? (Math.round(n.meta.alrt) + '/' + Math.round(n.meta.ufboot)) : Math.round(sv);
            if (note || n.meta.label) el('circle', { cx: x, cy: y, r: 4, fill: '#eda100', stroke: '#222', 'stroke-width': 1 }, g);
          }
          el('circle', { cx: x, cy: y, r: 6, fill: 'transparent' }, g); // hit target
          g.addEventListener('mousemove', ev => showTip(ev, `<b>${esc(n.meta.title || 'clade')}</b> · ${tipCount(n)} taxa` + (st ? `<br>${esc(st)}` : '') + (note ? `<br><i>${esc(note)}</i>` : '') + '<br><span class="mtree-hint">click to collapse</span>'));
          g.addEventListener('mouseleave', hideTip);
          g.addEventListener('click', () => { n.collapsed = true; draw(); });
        } else if (n.children.length && n.collapsed) {
          const cnt = tipCount(n), col = colorOf(n) || '#999', gname = groupOf(n) || 'mixed clade';
          const g = el('g', { class: 'mt-clade' }, gTips);
          const w = Math.min(60, 10 + Math.log2(cnt) * 8);
          el('path', { d: `M${x},${y}L${x + w},${y - th * 0.45}L${x + w},${y + th * 0.45}Z`, fill: col, 'fill-opacity': 0.55, stroke: col }, g);
          const label = el('text', { x: x + w + 5, y: y + fs / 3, fill: 'currentColor' }, g);
          label.textContent = `${gname} (${cnt})`;
          const st = supportText(n.meta);
          g.addEventListener('mousemove', ev => showTip(ev, `<b>${esc(gname)}</b> — ${cnt} taxa collapsed` + (st ? `<br>${esc(st)}` : '') + '<br><span class="mtree-hint">click to expand</span>'));
          g.addEventListener('mouseleave', hideTip);
          g.addEventListener('click', () => { n.collapsed = false; draw(); });
        } else {
          const name = displayName(n), col = n.meta.color || '#666';
          const hit = q && (name.toLowerCase().includes(q) || (n.meta.group || '').toLowerCase().includes(q) || n.label.toLowerCase().includes(q));
          if (hit) hits++;
          const g = el('g', { class: 'mt-tip' + (hit ? ' hit' : '') }, gTips);
          el('circle', { cx: x + 4, cy: y, r: 3.5, fill: col }, g);
          const t = el('text', { x: x + 10, y: y + fs / 3, fill: 'currentColor', 'font-style': n.meta.kind === 'mito' ? 'italic' : 'normal' }, g);
          t.textContent = name.length > 38 ? name.slice(0, 36) + '…' : name;
          const m = n.meta;
          const rows = [`<b>${esc(name)}</b>`];
          if (m.group) rows.push(esc(m.group) + (m.supergroup ? ' · ' + esc(m.supergroup) : ''));
          if (m.accession && m.accession !== name) rows.push('accession ' + esc(m.accession));
          if (m.taxid) rows.push('taxid ' + esc(m.taxid));
          if (m.note) rows.push('<i>' + esc(m.note) + '</i>');
          if (n.length !== null && phylogram) rows.push('branch length ' + Number(n.length).toPrecision(3));
          g.addEventListener('mousemove', ev => showTip(ev, rows.join('<br>')));
          g.addEventListener('mouseleave', hideTip);
          if (m.accession && /^(NC_|NW_|GCF_)/.test(m.accession)) {
            g.style.cursor = 'pointer';
            g.addEventListener('click', () => window.open(m.accession.startsWith('GCF_') ? 'https://www.ncbi.nlm.nih.gov/datasets/genome/' + m.accession : 'https://www.ncbi.nlm.nih.gov/nuccore/' + m.accession, '_blank'));
          }
        }
      })(root);
      // scale bar
      if (phylogram) {
        const nice = niceStep(L.maxD / 5), bw = nice / L.maxD * treeW;
        const gy = H - 12;
        el('path', { d: `M${padL},${gy}h${bw}`, stroke: 'currentColor', 'stroke-width': 1.2, fill: 'none' }, svg);
        el('text', { x: padL + bw / 2, y: gy - 4, 'text-anchor': 'middle', 'font-size': fs - 2, fill: 'currentColor' }, svg).textContent = nice + ' subst./site';
      }
      info.textContent = `${L.rows} rows · ${tipCount(root)} taxa` + (q ? ` · ${hits} match` : '');
      svg.style.color = getComputedStyle(container).color;
    }
    function niceStep(x) { const p = Math.pow(10, Math.floor(Math.log10(x))); const f = x / p; return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p; }
    function showTip(ev, htmlStr) {
      tooltip.innerHTML = htmlStr; tooltip.style.display = 'block';
      const r = container.getBoundingClientRect();
      let x = ev.clientX - r.left + 14, y = ev.clientY - r.top + 14;
      if (x + 260 > r.width) x = ev.clientX - r.left - 270;
      tooltip.style.left = x + 'px'; tooltip.style.top = y + 'px';
    }
    function hideTip() { tooltip.style.display = 'none'; }

    search.addEventListener('input', () => { query = search.value.trim(); draw(); });
    btnLen.addEventListener('click', () => { phylogram = !phylogram; btnLen.textContent = phylogram ? 'cladogram' : 'branch lengths'; draw(); });
    btnExpand.addEventListener('click', () => { (function w(n) { n.collapsed = false; n.children.forEach(w); })(root); draw(); });
    btnCollapse.addEventListener('click', () => {
      // collapse maximal single-group clades with ≥ 3 tips
      (function w(n) { if (!n.children.length) return; if (groupOf(n) && tipCount(n) >= 3 && n.parent) { n.collapsed = true; return; } n.children.forEach(w); })(root); draw();
    });
    btnSvg.addEventListener('click', () => {
      const svg = wrap.querySelector('svg'); if (!svg) return;
      const clone = svg.cloneNode(true); clone.setAttribute('style', 'color:#000;background:#fff');
      const blob = new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + clone.outerHTML], { type: 'image/svg+xml' });
      const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = (opts.filename || 'tree') + '.svg'; a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 1000);
    });
    if (opts.collapseGroups) btnCollapse.click(); else draw();
    let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(draw, 150); });
    return { redraw: draw, root };
  }

  async function load(container, url, opts) {
    try {
      const r = await fetch(url); if (!r.ok) throw new Error(r.status + ' ' + r.statusText);
      const data = await r.json();
      return render(container, data, opts);
    } catch (e) {
      container.innerHTML = '<p class="mtree-error">Could not load tree (' + esc(url) + '): ' + esc(e.message) + '</p>';
    }
  }

  global.MTree = { render, load, parseNewick };
})(window);
