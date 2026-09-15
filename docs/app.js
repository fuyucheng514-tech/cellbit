(() => {
  const pageOrder = ["home", "install", "usage", "input", "workflow", "stage3a", "stage3b", "advanced", "output", "api", "faq"];
  const pageNodes = [...document.querySelectorAll(".doc-page")];
  const pages = new Map(pageNodes.map((node) => [node.dataset.page, node]));
  const navLinks = [...document.querySelectorAll(".nav-link")];
  const sidebar = document.querySelector("#sidebar");
  const menuToggle = document.querySelector("#menuToggle");
  const searchInput = document.querySelector("#searchInput");
  const searchResults = document.querySelector("#searchResults");
  const searchEmpty = document.querySelector("#searchEmpty");
  const breadcrumb = document.querySelector("#breadcrumbCurrent");
  const tocList = document.querySelector("#tocList");
  const prevPage = document.querySelector("#prevPage");
  const nextPage = document.querySelector("#nextPage");

  function currentId() {
    const id = window.location.hash.slice(1).split("#")[0];
    return pages.has(id) ? id : "home";
  }

  function buildToc(page) {
    tocList.replaceChildren();
    page.querySelectorAll("h2[id], h3[id]").forEach((heading) => {
      const link = document.createElement("a");
      link.href = `#${page.dataset.page}#${heading.id}`;
      link.textContent = heading.textContent;
      if (heading.tagName === "H3") link.className = "toc-sub";
      link.addEventListener("click", (event) => {
        event.preventDefault();
        document.getElementById(heading.id)?.scrollIntoView({ block: "start" });
        history.replaceState(null, "", `#${page.dataset.page}#${heading.id}`);
      });
      tocList.appendChild(link);
    });
  }

  function setPage(id) {
    const safeId = pages.has(id) ? id : "home";
    pageNodes.forEach((page) => page.classList.toggle("active", page.dataset.page === safeId));
    navLinks.forEach((link) => link.classList.toggle("active", link.dataset.page === safeId));
    const page = pages.get(safeId);
    breadcrumb.textContent = page.dataset.title;
    buildToc(page);
    const index = pageOrder.indexOf(safeId);
    const previous = pageOrder[Math.max(0, index - 1)];
    const next = pageOrder[Math.min(pageOrder.length - 1, index + 1)];
    prevPage.href = `#${previous}`;
    prevPage.querySelector("span").textContent = index === 0 ? "Home" : "← Previous";
    prevPage.querySelector("strong").textContent = pages.get(previous).dataset.title;
    nextPage.href = `#${next}`;
    nextPage.querySelector("span").textContent = index === pageOrder.length - 1 ? "End" : "Next →";
    nextPage.querySelector("strong").textContent = pages.get(next).dataset.title;
    sidebar.classList.remove("open");
    menuToggle.setAttribute("aria-expanded", "false");
    document.title = `${page.dataset.title} · Microsags`;
    window.scrollTo({ top: 0, behavior: "instant" });
  }

  function closeSearch() {
    searchResults.hidden = true;
    searchResults.replaceChildren();
  }

  function search(query) {
    const needle = query.trim().toLowerCase();
    if (!needle) { closeSearch(); searchEmpty.hidden = true; return; }
    const matches = pageNodes.filter((page) => page.textContent.toLowerCase().includes(needle)).slice(0, 8);
    searchResults.replaceChildren();
    matches.forEach((page) => {
      const link = document.createElement("a");
      link.className = "search-result";
      link.href = `#${page.dataset.page}`;
      link.textContent = page.dataset.title;
      link.addEventListener("click", closeSearch);
      searchResults.appendChild(link);
    });
    searchResults.hidden = matches.length === 0;
    searchEmpty.hidden = matches.length !== 0;
  }

  navLinks.forEach((link) => link.addEventListener("click", () => { searchInput.value = ""; closeSearch(); }));
  window.addEventListener("hashchange", () => setPage(currentId()));
  searchInput.addEventListener("input", () => search(searchInput.value));
  document.addEventListener("keydown", (event) => {
    if (event.key === "/" && document.activeElement !== searchInput) { event.preventDefault(); searchInput.focus(); }
    if (event.key === "Escape") { closeSearch(); searchInput.blur(); }
  });
  document.addEventListener("click", (event) => { if (!event.target.closest(".search-wrap")) closeSearch(); });
  menuToggle.addEventListener("click", () => {
    const open = sidebar.classList.toggle("open");
    menuToggle.setAttribute("aria-expanded", String(open));
  });
  setPage(currentId());
})();
