# Beyond the Browser: demo materials

Materials from the AI session "Beyond the Browser: AI tools for better, faster research"
(Robert Hadfield, CRC 1211 IRTG Training Week, Bad Honnef, 8 October 2026).

Everything you need to see what the AI agent was given, what it produced, and to try it yourself.

## Links

- **Slides:** https://docs.google.com/presentation/d/13OPwpVx8OkNPJcNG9Z_TDxvvmtjVFq2nNoRRe84Cz2I/edit?usp=sharing
- **AI questionnaire results** (anonymous responses): https://docs.google.com/spreadsheets/d/1Zsx3o8ui1wtvBTezVloTgpM1qOYVOvnXcGhHdTY20tE/edit?usp=sharing
- **Science motivator** (the tool we built in the group work during the training week):
  https://radrob2.github.io/science-motivator/ (source: https://github.com/radrob2/science-motivator)
- **Free GitHub Copilot via GitHub Education** (setup steps): https://docs.google.com/document/d/1C_yjZ0_iAzP0AMjp3RJBskXo70ZPwiTVnGwONo0vEQU/edit?usp=sharing
- **Claude Team plan for scientists** (Anthropic): https://claude.com/programs/team-plan-for-scientists
  A PI (or equivalent) at a university or nonprofit research institute applies and verifies, then adds up to
  25 lab members. Standard seats are free (normally $20/month), premium seats with 5x the usage are
  $15/month (normally $100/month). Promotional pricing runs 12 months. PI verification form:
  https://claude.ai/labs-verification/attestation

## What's here

```
prompts/
  build_prompt.md          the prompt pasted into the agent to build the pipeline
  audit_prompt.md          the prompt pasted into a second, fresh agent session to audit the result
starting_project/          the project folder exactly as the agent first saw it
  docs/conventions.md      my project conventions, which the build prompt tells the agent to follow
  data/                    the demo input data (see "About the data" below)
  utils/normalize_fasta.py a small FASTA helper
  .claude/settings.json    pre-approves file edits so the agent does not stop to ask on every file
results/
  live_run_2026-10-08/     the live run during the session: config.sh, environment.yml, readme.md, src/,
                           share/ (tree, tables, report.md, summary.html), logs/, work/
  rehearsal_2026-10-02/    an earlier full rehearsal, plus TWO independent audits of it
                           (audit1/audit.md, audit2/audit.md, audit2/audit.html)
```

Open `results/live_run_2026-10-08/share/summary.html` in a browser to see the page projected in the session,
and `report.md` for the full write-up. The audit from the live run was not saved, but the exact prompt used is
`prompts/audit_prompt.md`. To see what an audit produces, open `results/rehearsal_2026-10-02/audit2/audit.html`.

The `config.sh`, `environment.yml`, `readme.md` and everything in `src/` were written by the agent. They are
not in `starting_project/`; producing them is part of the task.

## Try it yourself

You need an agent that can read and write files and run commands in a terminal (for example GitHub Copilot
CLI or Claude Code), conda or mamba, and an internet connection (the pipeline downloads from NCBI).

1. Copy `starting_project/` somewhere new and `cd` into it. Work on the copy, so you can start over.
2. Start your agent in that folder and paste in the contents of `prompts/build_prompt.md`.
3. When it finishes, open a **new, separate** agent session in the same folder and paste in
   `prompts/audit_prompt.md`. The auditor should not share memory with the builder.
4. Compare your results with `results/`. Your run will differ in the details; that is part of the point.

Notes:

- The build prompt asks for subagents and parallel work. Agents that cannot do this will just run steps in
  sequence.
- `docs/conventions.md` and the `data/.claude-sealed` marker refer to a sealed-data policy that is enforced
  on my machine (the agent may run tools on sealed data but not read it directly). That setup is not
  included here. With invented data it does not matter, but think about it before using real data.
- The rehearsal run used an earlier version of the build prompt without the five-minute compute cap. It
  pulled from both NCBI and BOLD, built the tree with IQ-TREE, and the whole build took about 32 minutes.
  The live run used NCBI only and FastTree, and its pipeline runs in a few minutes. The rehearsal
  folder's layout is flattened (no `share/` folder).

## About the data

**The input data is invented for this demo.** The sample names and sequences are not real Atacama
samples. Do not use any of the biological conclusions in `results/`.

<details>
<summary>Spoilers: what a good run should catch</summary>

- The ten "genome-extracted" sequences were made by lightly mutating a real public *Hyalella azteca*
  (North American) sequence. A careful run should flag them as implausible for the Atacama and suggest
  contamination or a sample mix-up, rather than report them as a finding.
- The build prompt says to keep the QC reference sequence out of the merged data, alignment and tree.
  In the rehearsal, the agent excluded it by name, but the same sequence came back in through the public
  download under its GenBank accession, and the report claimed it was absent. Both independent audits
  caught this by comparing sequence content. This is why the audit prompt asks to check by content, not
  only by name.

</details>

## License

Everything in this repository is dedicated to the public domain under [CC0 1.0](LICENSE). Use, copy, modify
and share it for any purpose, no attribution required. Public sequences downloaded from NCBI (in `work/`)
remain subject to NCBI's own terms.
