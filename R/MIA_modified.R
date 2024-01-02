MIA_modified=function(
    region_specific,               #cluster列是因子变量
    celltype_specific,             #cluster列是因子变量
    N,                             #所有基因的数量
    color_region,                  #配色，命名的字符串向量
    pvalue_log10_neg=20,           #p值相关的阈值
    axis.text.y.left.size=14,      #主panel纵轴的文本大小
    axis.ticks.y.left.length=0.2,  #主panel纵轴的刻度线长度
    legend.title.size=14,          #主panel的legend的title大小
    legend.ncol_uppanel = ceiling(length(unique(region_specific$region)) / 5), #上层panel的legend分为几列展示
    legend.text.size_uppanel=14,   #上层panel的legend的文本大小
    axis.text.y.left_uppanel=14    #上层panel的纵轴的文本大小
){
  miares=data.frame()
  for (ri in unique(as.character(region_specific$region))) {
    smalldf1=region_specific[region_specific$region == ri,]
    n=length(smalldf1$gene)
    term=c()
    pvalue=c()
    
    for (ci in unique(as.character(celltype_specific$celltype))) {
      one.celltype=celltype_specific[ celltype_specific$celltype %in% ci ,]
      M=length(one.celltype$gene)
      k=sum(smalldf1$gene %in% one.celltype$gene)
      one.pvalue=phyper(k-1,M, N-M, n, lower.tail=FALSE)
      term=append(term,ci)
      pvalue=append(pvalue,one.pvalue)
    }
    
    one.miares=data.frame(region=ri,term=term,pvalue=pvalue)
    one.miares=one.miares%>%arrange(pvalue)
    miares = miares%>%rbind(one.miares)
  }
  
  #添加几列信息
  miares$Enrichment= -log10(miares$pvalue)
  miares$Depletion= -log10(1-miares$pvalue)
  miares$final_value= ifelse(miares$Enrichment > miares$Depletion,miares$Enrichment,miares$Depletion)
  miares$final_class= ifelse(miares$Enrichment > miares$Depletion,"Enrichment","Depletion")
  miares$final_value[miares$final_value >= pvalue_log10_neg] = pvalue_log10_neg
  
  #画图部分
  labeldf1=as.data.frame(table(region_specific$region))
  colnames(labeldf1)=c("region","gene_num")
  #labeldf1$region=as.character(labeldf1$region)
  labeldf1$region_gene_num=paste0(as.character(labeldf1$region)," (",labeldf1$gene_num," genes)")
  labeldf1=labeldf1%>%arrange(region)
  
  labeldf2=as.data.frame(table(celltype_specific$celltype))
  colnames(labeldf2)=c("celltype","gene_num")
  #labeldf2$celltype=as.character(labeldf2$celltype)
  labeldf2$celltype_gene_num=paste0(as.character(labeldf2$celltype)," (",labeldf2$gene_num,")")
  labeldf2=labeldf2%>%arrange(celltype)
  
  miares$region=factor(miares$region,levels = levels(region_specific$region))
  miares$term=factor(miares$term,levels = levels(celltype_specific$celltype))
  miares.Enrichment=miares%>%filter(final_class == "Enrichment")
  miares.Depletion=miares%>%filter(final_class == "Depletion")
  
  library(ggnewscale)
  pb=ggplot()+
    geom_tile(data = miares.Enrichment,mapping = aes(x=region,y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Enrichment",colours = brewer.pal(9, "Reds"))+
    new_scale_fill() +
    geom_tile(data = miares.Depletion,mapping = aes(x=region,y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Depletion",colours = brewer.pal(9, "Blues"))+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0),breaks=labeldf2$celltype,labels=labeldf2$celltype_gene_num)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y.left = element_text(size = axis.text.y.left.size,color = "black"),
      axis.ticks.length.y.left = unit(axis.ticks.y.left.length,"cm"),
      axis.text.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = legend.title.size,vjust = 1)
    )
  
  updf=as.data.frame(levels(miares$region))
  colnames(updf)="region"
  updf$region=factor(updf$region,levels = updf$region)
  
  pa=ggplot(data = updf,aes(x=region,y=0,fill=region))+
    geom_tile()+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = color_region[as.character(updf$region)],breaks = labeldf1$region,labels=labeldf1$region_gene_num)+
    scale_y_continuous(expand = c(0,0),breaks = 0,labels = "Cell types (no. of genes)")+
    theme_void()+
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = legend.text.size_uppanel),
      legend.direction = "vertical",
      
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.text.y.left = element_text(size = axis.text.y.left_uppanel,color = "black")
    )+
    guides(fill = guide_legend(override.aes = list(size=10),ncol = legend.ncol_uppanel))
  
  library(patchwork)
  celltype_num=length(unique(celltype_specific$celltype))
  plot.f = pa / pb + plot_layout(heights = c(1,celltype_num))
  
  return(list(res.df = miares,res.plot = plot.f))
}
