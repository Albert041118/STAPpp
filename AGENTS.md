# AGENTS.md

## Role

You are working as a coding assistant for this repository.

Your priority is to complete the requested task with the smallest safe change.
Do not explore, refactor, or rewrite unrelated parts of the project unless explicitly asked.

## Working Rules

- First understand the specific task and the likely affected files.
- Prefer reading only the files directly related to the task.
- Do not scan the whole repository unless necessary.
- If more files are needed, briefly explain why before expanding the scope.
- Make minimal diffs.
- Do not reformat unrelated code.
- Do not rename files, move files, or change architecture unless required.
- Do not modify generated files, build outputs, or dependency lockfiles unless necessary.
- Do not add new dependencies unless the task cannot be solved reasonably without them.
- Preserve the existing code style and naming conventions.

## Before Editing

Before changing code, provide a short plan:

1. Suspected cause or target area
2. Files you plan to inspect or edit
3. Minimal change strategy
4. Relevant test or verification command

Do not start broad refactoring.

## During Editing

- Change only the files needed for the task.
- Keep the implementation simple.
- Avoid speculative improvements.
- Avoid touching unrelated warnings, lint issues, or TODOs.
- If you find a larger issue, report it separately instead of fixing it immediately.

## Testing

Run only the most relevant tests first.

Preferred order:

1. Unit test for the changed module
2. Related integration test
3. Type check or lint for affected files
4. Full test suite only if explicitly requested or clearly necessary

If a test fails because of unrelated existing issues, report that clearly and do not try to fix unrelated failures.

## Output Format

Keep responses concise.

After finishing, summarize:

- What changed
- Which files changed
- How it was verified
- Any remaining risks or follow-up suggestions

Avoid long explanations unless asked.
