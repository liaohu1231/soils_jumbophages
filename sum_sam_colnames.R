mean_same_col=function(x){
  t=table(colnames(x))
  nd<-matrix(0,length(t),ncol=length(rownames(x))) #得到一个矩阵,元素都为0,行数为t的个数，列数为d的变量个数减1
  for(i in 1:length(t)){
    index<-which(colnames(x) %in% names(t)[i])    #此行为匹配,把data的列名和t的第i个变量名进行匹配，并返回匹配到的下表值    
    if(length(index)==1) {
    nd[i,]<-x[1:length(rownames(x)),index]}   #当只匹配到一个值时，就不求总值了，直接存储就行了
    else 
    {nd[i,]<-rowSums(x[1:length(rownames(x)),index])}
    #rownames(nd)=rownames(data)
  }
  t=as.matrix(t)
  rownames(nd)=rownames(as.matrix(t))
  colnames(nd)=rownames(x)
  return(nd)
}
