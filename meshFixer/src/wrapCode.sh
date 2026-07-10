#!/bin/bash

OUTPUT="codebase_context.md"

# 1. Initialize file and generate the directory tree
echo "# Project Structure" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo '```text' >> "$OUTPUT"
# -a includes hidden files. -I ignores common noisy directories and the output file.
tree -a -I ".git|node_modules|venv|__pycache__|$OUTPUT" >> "$OUTPUT"
echo '```' >> "$OUTPUT"
echo -e "\n---\n" >> "$OUTPUT"

echo "# File Contents" >> "$OUTPUT"

# 2. Iterate through files and append their contents
find . -type f \
    -not -path '*/\.git/*' \
    -not -path '*/node_modules/*' \
    -not -path '*/venv/*' \
    -not -path '*/__pycache__/*' \
    -not -name "$OUTPUT" \
    | sort | while read -r file; do
    
    # Extract just the filename and the clean relative path
    filename=$(basename "$file")
    filepath="${file#./}"
    
    # Append the formatted data to the Markdown file
    echo "### File: $filename" >> "$OUTPUT"
    echo "**Location:** \`$filepath\`" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    echo '```' >> "$OUTPUT"
    cat "$file" >> "$OUTPUT"
    echo '```' >> "$OUTPUT"
    echo -e "\n---\n" >> "$OUTPUT"
done

echo "Done! Your context file has been generated: $OUTPUT"