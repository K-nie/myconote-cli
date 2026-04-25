// myconote compare report — vanilla JS, no dependencies.
// Two responsibilities: client-side table sort + ortholog filter.
// Keep this file under 500 lines per the v0.7.0 spec.

(function () {
  "use strict";

  // ── Sortable tables ─────────────────────────────────────────────────────
  function attachSort(table) {
    var headers = table.querySelectorAll("th[data-sortable='1']");
    headers.forEach(function (th, colIdx) {
      th.addEventListener("click", function () {
        var asc = !th.classList.contains("sorted-asc");
        // Clear sort markers on every header
        Array.prototype.forEach.call(headers, function (h) {
          h.classList.remove("sorted-asc", "sorted-desc");
        });
        th.classList.add(asc ? "sorted-asc" : "sorted-desc");
        sortRows(table, colIdx, asc);
      });
    });
  }

  function sortRows(table, colIdx, asc) {
    var tbody = table.tBodies[0];
    if (!tbody) return;
    var rows = Array.prototype.slice.call(tbody.rows);
    rows.sort(function (a, b) {
      var av = a.cells[colIdx].getAttribute("data-sort") ||
               a.cells[colIdx].textContent;
      var bv = b.cells[colIdx].getAttribute("data-sort") ||
               b.cells[colIdx].textContent;
      var an = parseFloat(av);
      var bn = parseFloat(bv);
      var cmp;
      if (!isNaN(an) && !isNaN(bn)) {
        cmp = an - bn;
      } else {
        cmp = av.localeCompare(bv);
      }
      return asc ? cmp : -cmp;
    });
    rows.forEach(function (r) { tbody.appendChild(r); });
  }

  // ── Ortholog table filter (text search) ─────────────────────────────────
  function attachFilter(input, tableId) {
    var table = document.getElementById(tableId);
    if (!table) return;
    input.addEventListener("input", function () {
      var q = input.value.toLowerCase();
      var rows = table.tBodies[0] ? table.tBodies[0].rows : [];
      Array.prototype.forEach.call(rows, function (row) {
        var hit = row.textContent.toLowerCase().indexOf(q) !== -1;
        row.style.display = (q === "" || hit) ? "" : "none";
      });
    });
  }

  // ── Bootstrap on DOMContentLoaded ───────────────────────────────────────
  document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll("table.sortable").forEach(attachSort);
    var search = document.getElementById("ortholog-search");
    if (search) {
      attachFilter(search, "ortholog-table");
    }
  });
})();
