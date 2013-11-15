if [ -z "$(svn status)" ]; then 
	LOCAL_CHANGES=""  
else 
	LOCAL_CHANGES="-Personal"
fi
echo $LOCAL_CHANGES
