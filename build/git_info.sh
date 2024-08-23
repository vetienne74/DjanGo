#!/usr/bin/env bash

git log | head -n3 > ./version.work
commit=$(grep -i commit ./version.work)
author=$(grep -i author ./version.work)
date=$(grep -i date ./version.work)
date_compile=$(date +%s)

#rm ./version.work

echo "const char DJANGO_GIT_COMMIT[] = \"$commit\" ;"
echo "const char DJANGO_GIT_AUTHOR[] = \"$author\" ;"
echo "const char DJANGO_GIT_DATE[] = \"$date\" ;"
echo "const time_t DJANGO_COMPILE_DATE = $date_compile ;"
