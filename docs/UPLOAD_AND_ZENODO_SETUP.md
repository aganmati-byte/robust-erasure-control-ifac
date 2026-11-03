# IFAC pre-publication repository — setup and access checklist (3 November 2026)

## Confirmed status

- GitHub account connected to ChatGPT: **aganmati-byte**.
- No separate IFAC GitHub repository existed among the repositories listed from the connected account during this preparation.
- A private repository with the name below is the recommended new home; the existing unrelated `adaptive-control-irregular-observations` repository must not be reused.
- A complete GitHub upload tree and a separate Zenodo data/software draft ZIP were prepared. **These files are local downloads, not yet uploaded to GitHub or Zenodo**.
- GitHub's current connector can inspect/edit *existing* repositories but cannot create a new repository or issue collaborator invitations from this chat.
- No Zenodo connector able to create an upload was available in this conversation. Zenodo draft creation requires your action in Zenodo's account interface (or a Work-mode cloud browser session).

## 1. Create private GitHub repository

1. Open https://github.com/new while signed into `aganmati-byte`.
2. Repository name: `robust-erasure-control-ifac`.
3. Description: `Private reproducibility research for robust mean-square stabilization under Bernoulli erasures (IFAC manuscript).`
4. Choose **Private**. Do not initialize with a README, .gitignore, or license (all are supplied in the upload tree).
5. Click **Create repository**.
6. Extract `IFAC_private_GitHub_repository_READY_2026-11-03.zip` on your computer. **Do not upload the ZIP as one file**: the repository should contain the actual folders/files.
7. On your new repository, click **Add file > Upload files**, drag the *extracted contents* into the uploader, then commit. GitHub supports folder uploads from a browser. Because this tree has fewer than 100 files, it fits GitHub's current browser limit.
8. Confirm that `README.md` is shown at the root, with `src/`, `code/`, `tests/`, `results/`, `figure4/` and `docs/`.
9. Settings > Collaborators > Add people. Invite Alexander using his GitHub username (not his academic email unless it is linked to his GitHub account), then wait for his acceptance. **Do not guess his GitHub login**.
10. Check the Actions tab for the CPU tests. The workflow is included in `.github/workflows/checks.yml`.

Alternative: create locally with GitHub Desktop, drag the extracted files into the repository folder, commit and select *Publish repository* with **Keep this code private** selected.

## 2. Prepare an **unpublished Zenodo draft** (do not publish)

1. Go to https://zenodo.org/ and sign in.
2. Select **New upload**. Resource type: **Software** (code and data supporting the IFAC manuscript).
3. Upload `IFAC_Zenodo_software_data_DRAFT_2026-11-03.zip`.
4. Use `zenodo/METADATA_DRAFT.json` from the GitHub package to fill the title, creators, description, and keywords. Correct both authors' affiliation wording with Alexander.
5. Choose a software/data license only after both authors agree. Do not invent grant details, ORCID identifiers, an article DOI or funding.
6. Click **Save draft**. **DO NOT CLICK PUBLISH.**
7. Click **Share** within the saved draft, give Alexander **Can edit** (if mutually agreed) or **Can preview draft** access; alternatively use a time-limited secret link with equivalent permission.
8. Reserve a DOI only when the record and version are stable and both authors agree; a reserved draft DOI is not yet a publicly published citable record.

## 3. Before any public deposit or IFAC submission

- Obtain Alexander's explicit review of the shared Overleaf PDF, full code, and mathematical/citation/novelty claims.
- Reconcile the **current Overleaf main.tex** with `manuscript/snapshot_2026-11-02_NOT_LIVE.tex`, which is only a dated copy.
- Confirm and document the exact Git commit/release identifier supporting the submitted figures and tables.
- Confirm authorship, corresponding-author role, affiliations, rights to release code, license, AI assistance, competing interests, funding and data statement.
- Check the journal's current Guide for Authors and supplement policy; do not assume Zenodo alone satisfies review-access requirements.
- Only after agreement decide whether to make the GitHub repo public and the Zenodo record published.

## Publishing strategy

GitHub = **living source-code collaboration**, Overleaf = **living manuscript**, Zenodo = **immutable, citable released code/data snapshot**. Draft Zenodo storage can remain unpublished while review continues. If preparing reviewer access, supply a journal-approved reviewer-accessible archive or time-limited secret link, not a private GitHub URL that reviewers cannot open.

## Source documentation

- GitHub official: https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository
- GitHub collaborator: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/repository-access-and-collaboration/inviting-collaborators-to-a-personal-repository
- Zenodo draft: https://help.zenodo.org/docs/deposit/create-new-upload/
- Zenodo draft sharing: https://help.zenodo.org/docs/share/user-sharing/
