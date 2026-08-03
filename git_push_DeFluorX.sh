#!/bin/bash
# =============================================================================
# DeFluorX - Git Push Helper
# Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
# Date   : 05 August 2026
# =============================================================================
# set -e is intentionally omitted: the script relies on explicit `$? -ne 0`
# checks after git commands. -u (unset-var guard) and pipefail are safe here.
set -uo pipefail

# ==========================================
# Terminal Colours & Formatting
# ==========================================
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
CYAN='\033[1;36m'
NC='\033[0m'

REPO="KU-MGB/DeFluorX"

FORCE=false
for arg in "$@"; do
    case "$arg" in
        -f|--force)
            FORCE=true
            ;;
    esac
done

echo -e "${BLUE}==========================================${NC}"
echo -e "${CYAN}  Git Push - ${REPO}${NC}"
echo -e "${BLUE}==========================================${NC}\n"

# ==========================================
# 0. Mode notice - give the user a chance to abort and switch to force
# ==========================================
if [ "$FORCE" = true ]; then
    echo -e "${RED}FORCE MODE: the online GitHub repository will be WIPED and replaced to exactly match your local directory.${NC}"
    read -p "Press ENTER to continue, or Ctrl+C to cancel... " _
else
    echo -e "${YELLOW}Normal push mode.\nTo completely wipe the online GitHub repository and duplicate your local directory, cancel now (Ctrl+C) and re-run with --force or -f.${NC}"
    read -p "Press ENTER to continue with a normal push, or Ctrl+C to cancel... " _
fi
echo

# ==========================================
# 1. Status & Branch Info
# ==========================================
echo -e "${BLUE}Checking Git status...${NC}"
git status -s

current_branch=$(git branch --show-current)
# Detached HEAD → --show-current is empty; a later `git push origin ""` fails fatally.
if [ -z "$current_branch" ]; then
    echo -e "\n${RED:-}ERROR: detached HEAD - no branch to push. Check out a branch first (e.g. 'git switch main').${NC:-}" >&2
    exit 1
fi
echo -e "\n${CYAN}Branch: ${NC}${current_branch}\n"

# ==========================================
# 2. Staging & Large File Check
# ==========================================
echo -e "${BLUE}Adding all tracked changes...${NC}"
git add -A

echo -e "${BLUE}Checking staged files for large blobs (>100MB)...${NC}"
# Only inspect files actually staged for commit. This respects .gitignore /
# .git/info/exclude for free (ignored paths are never staged) and avoids
# walking huge ignored trees like Boltz-2_*/, __pycache__/, Older_Codes_Zipped/.
# The size that matters is the STAGED BLOB's, not the working tree's: a file can be staged large
# and then shrunk (or deleted) on disk, and GitHub would still reject the pushed object. Ask the
# object database for the size of what is actually staged.
# -z / read -d '': paths are NUL-separated, so a filename containing a space or a newline is read
# as one path instead of being split into fragments that then fail every lookup silently.
large_files=$(git diff --cached -z --name-only --diff-filter=ACM | while IFS= read -r -d '' f; do
    sha=$(git ls-files -s -- "$f" | awk '{print $2}')
    [ -n "$sha" ] || continue
    sz=$(git cat-file -s "$sha" 2>/dev/null || echo 0)
    [ "$sz" -gt 104857600 ] && printf '%s (%s bytes, staged)\n' "$f" "$sz"
done)

if [ -n "$large_files" ]; then
    echo -e "${RED}ERROR: Large files (>100MB) are staged and would be rejected by GitHub:${NC}"
    echo -e "${YELLOW}$large_files${NC}"
    echo -e "${RED}Aborting before commit. Unstage/remove them (or use Git LFS), then re-run.${NC}\n"
    exit 1
else
    echo -e "${GREEN}No large files detected.${NC}\n"
fi

# ==========================================
# 3. Check if there is anything to commit
# ==========================================
if git diff --quiet && git diff --cached --quiet; then
    echo -e "${GREEN}Nothing to commit. Working tree clean.${NC}"
    exit 0
fi

# ==========================================
# 4. Python Syntax Check on Staged Files
# ==========================================
if git diff --cached --name-only --diff-filter=ACM | grep -q '\.py$'; then
    echo -e "${BLUE}Running syntax check on staged Python files...${NC}"
    git diff --cached -z --name-only --diff-filter=ACM | while IFS= read -r -d '' file; do
        case "$file" in *.py) ;; *) continue ;; esac
        # Validate the STAGED blob (git show :file), not the working-tree copy - otherwise a
        # broken file that was staged then fixed on disk would wrongly pass.
        if ! git show ":$file" | python3 -c "import ast,sys; ast.parse(sys.stdin.read())" 2>/dev/null; then
            echo -e "\n${RED}Syntax error in staged: $file - fix before committing. Aborting.${NC}"
            exit 1
        fi
    done || exit 1
    echo -e "${GREEN}Python syntax check passed.${NC}\n"
fi

# ==========================================
# 4b. Bash Syntax Check on Staged Shell Files
# ==========================================
if git diff --cached --name-only --diff-filter=ACM | grep -q '\.sh$'; then
    echo -e "${BLUE}Running syntax check on staged shell scripts...${NC}"
    git diff --cached -z --name-only --diff-filter=ACM | while IFS= read -r -d '' file; do
        case "$file" in *.sh) ;; *) continue ;; esac
        # Validate the STAGED blob, not the working-tree copy: a script that was staged broken and
        # then fixed on disk would otherwise pass the check and be pushed broken. Matches the
        # Python check above.
        if ! git show ":$file" | bash -n /dev/stdin; then
            echo -e "\n${RED}Syntax error in staged: $file - fix before committing. Aborting.${NC}"
            exit 1
        fi
    done || exit 1
    echo -e "${GREEN}Shell syntax check passed.${NC}\n"
fi

# ==========================================
# 5. Conventional Commit Prompt
# ==========================================
echo -e "${CYAN}Select commit type:${NC}"
echo "1) revamp (complete code revamp)"
echo "2) feat   (new script or feature)"
echo "3) fix    (bug fix or correction)"
echo "4) docs   (notes, README, text files)"
echo "5) data   (structures, PDBs, datasets)"
echo "6) chore  (maintenance, moving files)"
read -p "> " type_choice

case $type_choice in
    1) c_type="revamp" ;;
    2) c_type="feat" ;;
    3) c_type="fix" ;;
    4) c_type="docs" ;;
    5) c_type="data" ;;
    6) c_type="chore" ;;
    *) c_type="update" ;;
esac

echo -e "\n${CYAN}Enter commit message:${NC}"
read -p "> " msg

if [ -z "$msg" ]; then
    echo -e "${RED}Commit message cannot be empty. Aborting.${NC}"
    exit 1
fi

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
full_msg="[$timestamp] $c_type: $msg"

# ==========================================
# 6. Branch Protection & Push
# ==========================================
if [ "$current_branch" == "main" ] || [ "$current_branch" == "master" ]; then
    echo -e "\n${YELLOW}You are pushing directly to '${current_branch}'.${NC}"
    read -p "Continue? (y/n): " confirm
    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
        echo -e "\n${RED}Push aborted.${NC}"
        exit 0
    fi
fi

echo -e "\n${BLUE}Committing...${NC}"
git commit -m "$full_msg" --quiet

if [ $? -ne 0 ]; then
    echo -e "${RED}Commit failed.${NC}"
    exit 1
fi

echo -e "\n${BLUE}Pushing to ${REPO} (branch: ${current_branch})...${NC}"
if [ "$FORCE" = true ]; then
    echo -e "${CYAN}Remote will be fast-forwarded to match your local branch (force-with-lease: aborts if the remote moved since your last fetch).${NC}"
    git push --force-with-lease --set-upstream origin "$current_branch"
else
    git push --set-upstream origin "$current_branch"
fi

if [ $? -ne 0 ]; then
    echo -e "${RED}Push failed. Check your connection or remote.${NC}"
    exit 1
fi

short_hash=$(git rev-parse --short HEAD)
echo -e "\n${GREEN}Done! Pushed to ${REPO}.${NC}"
echo -e "${CYAN}Link: https://github.com/${REPO}/commit/${short_hash}${NC}"
echo -e "${BLUE}==========================================${NC}"
