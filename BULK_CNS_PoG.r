library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggiraph)
library(qusage)
library(RColorBrewer)
library(ggpubr)
library(shinythemes)
library(ggplotify)
library(edgeR)
library(org.Hs.eg.db)
library(ggrepel)
library(nloptr)
library(BiocManager)

exp = readRDS('bulk_analysis.rds')
load("res.RData")

dat=read.csv('merged_counts.csv',sep=';', header = TRUE, row.names=1)
meta=read.csv('meta_data_3batch.csv', header = T,sep = ';', row.names = 1)
sign_cols=c(DN='deepskyblue',NS='grey80',UP='tomato')
res=rbind(resA,resB,resC,resD,resE,resF,resG,resH,resI,resJ)
res= na.omit(res)
#dat=exp$norm.corrected
#meta=exp$meta
conv=exp$conv
#rownames(dat)=conv[rownames(dat),2]
#names(meta)[names(meta) == "type"] <- "type"
meta$sex = 'M'
meta[which(meta$sample=="AB203"),]$sex = 'F'
meta[which(meta$sample=="AB253"),]$sex = 'F'
meta$condition = 'Normal'
meta[which(meta$type %in% c("WML","CL","PVML")),'condition']="Lesion"
meta[grepl('^NA',meta$type),'condition']='Normal_Appearing'

d0 <- DGEList(dat[-1])
d0 <- calcNormFactors(d0)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
mm <- model.matrix(~ 0 + type+batch+sample, data=meta[colnames(d),])
y <- voom(d[,rownames(meta)], mm, plot = F)
dat <- removeBatchEffect(y$E,batch = meta$batch, batch2 = meta$sample)
rownames(dat)=conv[rownames(dat),2]
meta$gene = dat['MCAM',rownames(meta)]

plot.geneDEA <- function(x){
  gene_id = rownames(subset(resA,gene_name==x))[1]
  meta$gene = y$E[gene_id,rownames(meta)]
  ggplot(meta,aes(type,gene))+geom_boxplot()+geom_point(aes(col=sample,shape=batch))+theme_bw()+ ggtitle(x)
}

plot.genePCA<-function(x,group='type',meta){
    meta$gene=dat[x,rownames(meta)]
    
    if(!is.numeric(meta[,group])){
        t=table(meta[,group])
        ggplot(meta[which(meta[,group]%in%names(t[t>2])),],aes(meta[which(meta[,group]%in%names(t[t>2])),group],gene,fill=meta[which(meta[,group]%in%names(t[t>2])),group]))+geom_boxplot(outlier.shape=NA)+guides(fill=F)+ggtitle(x)+theme_bw()+geom_point(aes(col=meta[which(meta[,group]%in%names(t[t>2])),'batch']),position=position_jitter(width=.1))+xlab('')+labs(col='sample')+theme(axis.text.x=element_text(angle=60,hjust=1))+ylab('')+scale_color_manual(values=c('gray10','gray60','gray80'))+stat_compare_means()
    }else{
        ggplot(meta,aes(meta$gene,meta[,group]))+geom_point()+theme_bw()+stat_cor()
    }
}

go = read.gmt('c5.bp.v7.0.symbols.gmt')
names(go)=gsub('_',' ',gsub('GO_','',names(go)))

ui <- navbarPage("Basic dashboard", theme = shinytheme("flatly"),
  tabPanel(
      'PCA',
      sidebarLayout(
             sidebarPanel(width=1),           
    mainPanel(

    fluidRow(
      box(
          width=12,
          column(width=4,radioButtons("geneset", 'Geneset', c('All','GO'), selected = 'All',inline=T)),
          column(width=6,selectizeInput('color_by','Color by',c('batch','type','tissue','condition','sample','sex','expression'),selected='batch')),
          conditionalPanel(
              "input.color_by=='expression'",
              column(width=4,selectizeInput('gene','Gene',rownames(dat),'MCAM'))
          ),
          conditionalPanel(
              "input.geneset=='HVG'",
              column(width=6,sliderInput('n_hvg','Nb HVG',min=10,max=nrow(dat),value=3000,step=5))
          ),
          conditionalPanel(
              "input.geneset=='GO'",
              column(width=6,selectizeInput('go_term','GO Term',names(go),selected='CELL ADHESION MOLECULE PRODUCTION'))
          )),
        box(width=12,girafeOutput('reduction')),
        conditionalPanel(
            "input.color_by=='expression'",
            box(
                width=12,
                selectizeInput('group_by','Group by',c('batch','type','tissue','condition','sample','sex','DIM_1','DIM_2'),selected='batch'),
                plotOutput('boxplot')
            )
        )
      
    )
  )
)),
tabPanel('DEA',
         inputPanel(selectInput('volcano','Select Comparison:',choices = unique(res$comp),selected = "PVML_vs_NAWM"),
                    numericInput("top", "Number of top genes:", 15,
                                 min = 5, max = 100),
                    selectInput('geneb','Select Top genes by:',choices = c("large fold changes","high statistical significance")), selectInput('GENE','Select gene:',as.character(unique(res$gene_name))),checkboxInput("topshow", "Show top genes", TRUE), 
                    checkboxInput("Vplot", "Multiple view (All VolcanoPlot)", FALSE)),
         fluidRow(
           column(7, 
                  girafeOutput("PlotVolcano")
           ),
           column(5, 
                  plotOutput("PlotGene")
           )
         ) 
))

server <- function(input, output) {
  set.seed(122)
  histdata <- rnorm(500)

  output$reduction <- renderGirafe({
    
    #num = dat
    if(input$geneset=='All'){
        pca = prcomp(t(dat))
        
    }else if(input$geneset=='GO'){
        set=go[[input$go_term]]
        num=dat[set[set%in%rownames(dat)],]
        pca = prcomp(t(num))
    }
    meta$DIM_1 =  pca$x[rownames(meta),1]
    meta$DIM_2 =  pca$x[rownames(meta),2]
    if(input$color_by=='expression'){
        meta$gene = dat[input$gene,rownames(meta)]
        color_by='gene'
    }else{
        color_by=input$color_by
    }
    pl = ggplot(meta,aes(DIM_1,DIM_2,fill=meta[,color_by]))+geom_point_interactive(aes(tooltip=id,data_id=id),shape=21,size=3)+theme_bw()
    if(input$color_by=='expression'){
        pl=pl+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))
    }
    girafe(ggobj=pl+labs(fill=input$color_by))
  })

  output$boxplot<-renderPlot({
    if(input$geneset=='All'){
        pca = prcomp(t(dat))
        
    }else if(input$geneset=='GO'){
        set=go[[input$go_term]]
        num=dat[set[set%in%rownames(dat)],]
        pca = prcomp(t(num))
    }
      meta$DIM_1 =  pca$x[rownames(meta),1]
      meta$DIM_2 =  pca$x[rownames(meta),2]
      pl=plot.genePCA(input$gene,input$group_by,meta)
      pl
  })

  output$PlotVolcano<-renderGirafe({
    if(input$geneb=="large fold changes"){
      resA=resA[order(-abs(resA$logFC)),]
      resB=resB[order(-abs(resB$logFC)),]
      resC=resC[order(-abs(resC$logFC)),]
      resD=resD[order(-abs(resD$logFC)),]
      resE=resE[order(-abs(resE$logFC)),]
      resF=resF[order(-abs(resF$logFC)),]
      resG=resG[order(-abs(resG$logFC)),]
      resH=resH[order(-abs(resH$logFC)),]
      resI=resI[order(-abs(resI$logFC)),]
      resJ=resJ[order(-abs(resJ$logFC)),]
      
    }else if(input$geneb=="high statistical significance"){
      resA=resA[order(abs(resA$P.Value)),]
      resB=resB[order(abs(resB$P.Value)),]
      resC=resC[order(abs(resC$P.Value)),]
      resD=resD[order(abs(resD$P.Value)),]
      resE=resE[order(abs(resE$P.Value)),]
      resF=resF[order(abs(resF$P.Value)),]
      resG=resG[order(abs(resG$P.Value)),]
      resH=resH[order(abs(resH$P.Value)),]
      resI=resI[order(abs(resI$P.Value)),]
      resJ=resJ[order(abs(resJ$P.Value)),]
    }
    resA$show = F
    resA[1:input$top,'show']=T
    resB$show = F
    resB[1:input$top,'show']=T
    resC$show = F
    resC[1:input$top,'show']=T
    resD$show = F
    resD[1:input$top,'show']=T
    resE$show = F
    resE[1:input$top,'show']=T
    resF$show = F
    resF[1:input$top,'show']=T
    resG$show = F
    resG[1:input$top,'show']=T
    resH$show = F
    resH[1:input$top,'show']=T
    resI$show = F
    resI[1:input$top,'show']=T
    resJ$show = F
    resJ[1:input$top,'show']=T
    
    res=rbind(resA,resB,resC,resD,resE,resF,resG,resH,resI,resJ)
    res= na.omit(res)
    res$sign = 'NS'
    res[which(res$adj.P.Val<0.05&res$logFC>1),'sign']='UP'
    res[which(res$adj.P.Val<0.05&res$logFC<(-1)),'sign']='DN'
    
    if(input$Vplot){
      pl = ggplot(res,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,data_id = gene_name,data_id = gene_name))+facet_wrap(~comp, ncol=4)+theme_bw()+scale_color_manual(values=sign_cols)
      if(input$topshow){
        pl = pl+geom_label_repel(data=subset(res,show),aes(label=gene_name),col='black',size=1.75,force=30,fill = alpha(c("white"),0.5))
      }
    }
    else{ 
      pl = ggplot(subset(res,res$comp==input$volcano),aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)
      if(input$topshow){
        pl = pl+geom_label_repel(data=subset(subset(res,res$comp==input$volcano),show),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))
      }
    }
    x = girafe(ggobj = pl,width_svg = 8, height_svg = 8)
    x = girafe_options(x,opts_selection(type="single"))
    return(x)
  })
  
  values = reactiveValues()
  values$show = 'CD24'
  
  observeEvent(input$GENE,{
    values$show = input$GENE
  })
  
  output$PlotGene <- renderPlot({plot.geneDEA(values$show)},height = 300, width = 450)
}

shinyApp(ui, server)