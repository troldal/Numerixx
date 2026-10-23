# reviews/

One file per review: `R<round>-<reviewer>.md` (for example `R1-numerics.md`), written by the reviewer agent, never
edited afterwards. Each finding has an id (`R1-numerics-07`), a location (chapter file and line, or figure, table,
listing), a severity (major, minor, editorial), the problem with its evidence, and a suggested change. The editor
writes `R<round>-decision.md`: the decision on every finding (accept, reject with the reason, or defer) and the
revision plan by chapter. The independent review by another model (OpenAI Codex) is `R<round>-codex.md`.
