#对x、y轴name做修改
p+theme(axis.title= element_text(size=15, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5))
#仅对y轴name做修改：
p+theme(axis.title.y= element_text(size=15, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5))
#仅对y轴的刻度做修改：
p+theme(axis.text.y= element_text(size=15, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5))
#对图片的title做修改
p+theme(title= element_text(size=15, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5))
#对legend的内容做修改
p+theme(legend.text= element_text(size=15, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5))