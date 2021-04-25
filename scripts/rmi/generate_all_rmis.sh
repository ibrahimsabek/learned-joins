rm -f $(dirname "$0")/../../include/rmi/all_rmis.h

# An adapted version of gen.sh in the SOSD benchmark 
echo "#pragma once" >> $(dirname "$0")/../../include/rmi/all_rmis.h
for header in $(ls $(dirname "$0")/../../include/rmi/ | grep "\\.h$" | grep -v data | grep -v all_rmis ); do
    echo "#include \"${header}\"" >> $(dirname "$0")/../../include/rmi/all_rmis.h
done
