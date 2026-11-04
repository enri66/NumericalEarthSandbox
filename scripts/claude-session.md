# Claude on antares over SSH

`claude-session` runs Claude Code inside tmux on antares, so it survives SSH drops and keeps working while you're away.

## Install

```bash
cd ~/Models/NumericalEarthSandbox && git pull origin main
mkdir -p ~/bin && ln -s ~/Models/NumericalEarthSandbox/scripts/claude-session ~/bin/claude-session
```

If `claude` isn't installed yet: `curl -fsSL https://claude.ai/install.sh | bash`, then run `claude` once to log in (it prints a URL to open on your laptop or phone, and you paste the code back).

## The loop

1. `ssh antares`
2. `claude-session` opens a tmux session in `~/Models/NumericalEarthSandbox` at the Claude prompt, continuing the last conversation there if there is one.
3. Type what you want, e.g. *"check how the mab 90-day run is doing: tail the driver log, count SPIKEs and NaNs"*. It runs on antares, so it can read logs, run `squeue`/`pgrep`, edit scripts, launch Julia and push to GitHub itself.
4. When it asks to run a command or edit a file, Enter approves (or choose "don't ask again" for trusted commands). Esc stops it mid-step.
5. Done for now: **Ctrl-b d** (detach). Claude keeps running in tmux; you can close the laptop or let SSH drop.
6. Next time: `ssh antares`, `claude-session`, and you're back on the same screen with the history intact.

Don't close the terminal or type `exit` to leave; use Ctrl-b d.

## Options

| Command | What it does |
|---|---|
| `claude-session` | attach, or start and continue the last conversation |
| `claude-session -n` | start a fresh conversation |
| `claude-session -r` | run `claude remote-control`, so the session also shows up in the Claude app on your phone or laptop |
| `claude-session -s NAME` | use tmux session `NAME` (one per project works well) |
| `claude-session -d DIR` | working directory (default: `$CLAUDE_SESSION_DIR` or `~/Models/NumericalEarthSandbox`) |

`-n` and `-r` only take effect when the tmux session doesn't exist yet. To switch a running one to `-r`: `/exit` in Claude, `exit` the shell that's left, then `claude-session -r`.

## Inside Claude

| Command | What it does |
|---|---|
| `/resume` | switch to an earlier conversation |
| `/clear` | start fresh on a new topic |
| `/exit` | quit Claude (you land on a shell in that tmux window; type `claude` to restart) |
