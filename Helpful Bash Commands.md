# Helpful Bash Commands Cheatsheet

A collection of commonly used **Bash / Linux shell commands**.

---

## File & Directory Management

```bash
ls            # List files
ls -la        # List all files (including hidden) with details
pwd           # Show current directory
cd /path      # Change directory
cd ..         # Go up one level
cd ~          # Go to home directory
mkdir dir     # Create directory
mkdir -p a/b  # Create nested directories
rmdir dir     # Remove empty directory
rm file       # Remove file
rm -r dir     # Remove directory recursively
rm -rf dir    # Force remove (⚠ dangerous)
cp src dst    # Copy file
cp -r a b     # Copy directory
mv src dst    # Move or rename
stat file     # File details
```

---

## File Viewing & Editing

```bash
cat file            # Print file contents
less file           # Paginated view (recommended)
head file           # First 10 lines
head -n 50 file     # First 50 lines
tail file           # Last 10 lines
tail -f logfile     # Follow file updates
nano file           # Simple terminal editor
vim file            # Advanced editor
```

---

## Searching & Finding

```bash
grep "text" file           # Search text in file
grep -r "text" dir         # Recursive search
find . -name "*.sh"        # Find by name
find . -type f             # Find files only
find . -type d             # Find directories
which cmd                  # Command location
whereis cmd                # Command binaries
```

---

## Compression & Archives

```bash
gzip file                  # Compress
unzip file.zip             # Extract zip
zip -r file.zip dir        # Zip directory
tar -cvf a.tar dir          # Create tar
tar -xvf a.tar              # Extract tar
tar -czvf a.tar.gz dir     # Create tar.gz
tar -xzvf a.tar.gz         # Extract tar.gz
```

---

## Permissions & Ownership

```bash
chmod +x file               # Make executable
chmod 755 file              # rwxr-xr-x
ls -l                       # View permissions
```

---

## Networking

```bash
wget url                     # Download file
```

---

## Processes & System Info

```bash
ps aux                       # List processes
top                          # Live process viewer
htop                         # Better top (if installed)
kill PID                     # Kill process
kill -9 PID                  # Force kill
free -h                      # Memory usage
df -h                        # Disk usage
```

---

## Input / Output Redirection

```bash
cmd > file        # Redirect output
cmd >> file       # Append output
cmd < file        # Input from file
cmd1 | cmd2       # Pipe output
```

---

## Environment Variables

```bash
printenv               # List env vars
echo $PATH             # Print variable
export VAR=value       # Set variable
unset VAR              # Remove variable
source file.sh         # Load variables
```

---

## Package Management

### Debian / Ubuntu

```bash
sudo apt update
sudo apt install pkg
sudo apt remove pkg
```


---


