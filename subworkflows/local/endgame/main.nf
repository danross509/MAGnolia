
    """
        for f in \$(find . -path ./work -prune -o -type l -print); do
            target=\$(readlink -e \$f)
            sed -i '' \$f
            rm \$target
        done
    """