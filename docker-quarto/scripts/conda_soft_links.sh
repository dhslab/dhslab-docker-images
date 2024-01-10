# this loop to avoid overiding any existing soft-linked binary
link_path=$1
for file in `ls $link_path` ; do 
    if [ -f /usr/local/bin/$file ] 
    then 
        :
    else 
        ln -s ${link_path}/$file /usr/local/bin/$file
    fi 
done