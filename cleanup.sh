# Check current disk usage
echo "🔍 Disk usage:"
df -h ~

# Clear user trash folders
echo "🗑 Clearing Trash..."
rm -rf ~/.local/share/Trash/* ~/.Trash/* 2>/dev/null

# Clear Jupyter runtime files
echo "🧹 Clearing Jupyter runtime..."
rm -rf ~/.local/share/jupyter/runtime/*

# Clear IPython caches
echo "🧹 Clearing IPython cache..."
rm -rf ~/.ipython/profile_default/db
rm -rf ~/.ipython/profile_default/security

# Clear temp files (only your user-owned)
echo "🧹 Clearing /tmp files owned by you..."
find /tmp -user "$USER" -type f -delete 2>/dev/null

echo "✅ Cleanup complete."
