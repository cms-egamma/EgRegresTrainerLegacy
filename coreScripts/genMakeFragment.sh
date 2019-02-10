package=$1

echo $package | awk '{printf "%s%s",toupper($1),"_LIBFILES\t="}'

for file in `ls packages/$package/src/*.cc`
do
  echo $file | awk -F "/" '{printf "%s%s ","$(PKG_OBJ_DIR)/",$(NF)}' | sed 's/.cc/.o/g'
done
#echo ""


for file in `ls packages/$package/dict/*_LinkDef.h`
do
  echo $file | awk -F "/" '{printf "%s%s ","$(PKG_OBJ_DIR)/",$(NF)}' | sed 's/_LinkDef.h/Dict.o/g'
done
echo ""