library("igraph")
library("dplyr")

load("mic_meta_cor_r.rda")
load("mic_meta_cor_p.rda")
load("node.rda")
load("differ.rda")
load("met_abun.rda")
load("mic_abun.rda")
load("met_diff.rda")
load("mic_diff.rda")


#--------------- Parameter Settings

# Threshold of correlation coefficient and p value
r_threshold <- 0.9
p_threshold <- 0.05


#--------------- Generate edge table
CorrDF <- function(cormat, pmat) {
  ut = matrix (TRUE, nrow = nrow (cormat), ncol = ncol (pmat))
  data.frame(
    from = rownames(cormat)[row(cormat)[ut]],
    to = colnames(cormat)[col(cormat)[ut]],
    r =(cormat)[ut],
    p.adjust = pmat[ut]
  )
}
mic_meta_cor_df <- CorrDF(mic_meta_cor_r,mic_meta_cor_p)
# Screening edge
mic_meta_cor_filtered <- mic_meta_cor_df[which(abs(mic_meta_cor_df$r) > r_threshold & mic_meta_cor_df$p.adjust < p_threshold),]


#--------------- Construct a correlation network
nodeattrib <- data.frame(Node = union(mic_meta_cor_filtered$from,mic_meta_cor_filtered$to))
nodeattrib <- merge(nodeattrib,node,all.x = T)
rownames(nodeattrib) <- nodeattrib$Node

# Draw a co-occurrence network diagram with igraph
co_net <- graph_from_data_frame(mic_meta_cor_filtered,direct = F, vertices = nodeattrib) 

# Nodes attribute table
net_nodes <- rownames(nodeattrib[nodeattrib$group %in% co_net,])




# I. Calculate degree score(DS)
{
  nodes_size = centr_degree(co_net)$res
  DS <- as.data.frame(cbind(names(V(co_net)),nodes_size))
  colnames(DS) <- c("Node","size")
  DS$size <- as.numeric(DS$size)
  DS <- DS[order(DS$size,decreasing = T),]
  DS$score <-  seq(nrow(DS),1,-1)
  DS <- merge(DS,node,all.x = T)

  DS_metabolites <- DS[which(DS$Type == "metabolite"),]
  DS_met_function <- DS[which(DS$Type == "metabolic_function"),]
  DS_microbes <- DS[which(DS$Type == "microbe"),]
  DS_mic_function <- DS[which(DS$Type == "microbial_function"),]
  
  # Convert the score to[1,3]
  DS_metabolites$score <-  seq(nrow(DS_metabolites),1,-1)
  DS_metabolites$score <- (DS_metabolites$score - 1)/(nrow(DS_metabolites) - 1) * 2 + 1
  DS_met_function$score <-  seq(nrow(DS_met_function),1,-1)
  DS_met_function$score <- (DS_met_function$score - 1)/(nrow(DS_met_function) - 1) * 2 + 1
  DS_microbes$score <-  seq(nrow(DS_microbes),1,-1)
  DS_microbes$score <- (DS_microbes$score - 1)/(nrow(DS_microbes) - 1) * 2 + 1
  DS_mic_function$score <-  seq(nrow(DS_mic_function),1,-1)
  DS_mic_function$score <- (DS_mic_function$score - 1)/(nrow(DS_mic_function) - 1) * 2 + 1
}


# II. Calculate edge score(ES)
{
  # If both nodes are differential nodes, ES is equal to 3 points, if only one of them is differential, ES is equal to 2, and if neither of them is differential, ES is equal to 1
  mic_meta_cor_filtered$ES <- 1
  dfii <- unique(differ$Diff)
  for (i in 1:nrow(mic_meta_cor_filtered)) {
    if((mic_meta_cor_filtered$from[i] %in% dfii) & (mic_meta_cor_filtered$to[i] %in% dfii)){
      mic_meta_cor_filtered$ES[i] <- 3
    }else if((mic_meta_cor_filtered$from[i] %in% dfii) | (mic_meta_cor_filtered$to[i] %in% dfii)){
      mic_meta_cor_filtered$ES[i] <- 2
    }else{
      mic_meta_cor_filtered$ES[i] <- 1
    }
  }

  from_es <- mic_meta_cor_filtered %>%
    group_by(from) %>%
    summarise("ES" = sum(ES))
  colnames(from_es) <- c("Node","ES")

  to_es <- mic_meta_cor_filtered %>%
    group_by(to) %>%
    summarise("ES" = sum(ES))
  colnames(to_es) <- c("Node","ES")
  
  mic_met_es <- merge(from_es,to_es,by = "Node",all = T)
  mic_met_es[is.na(mic_met_es)] <- 0
  mic_met_es$ES <- mic_met_es$ES.x + mic_met_es$ES.y
  mic_met_es <- mic_met_es[,-(2:3)]
  ES <- mic_met_es
  ES <- ES[order(ES$ES,decreasing = T),]
  ES$score <-  seq(nrow(ES),1,-1)
  ES <- as.data.frame(ES)
  ES <- merge(ES,node,all.x = T)
  
  ES_metabolites <- ES[which(ES$Type == "metabolite"),]
  ES_met_function <- ES[which(ES$Type == "metabolic_function"),]
  ES_microbes <- ES[which(ES$Type == "microbe"),]
  ES_mic_function <- ES[which(ES$Type == "microbial_function"),]
  
  # Convert the score to[1,3]
  ES_metabolites$score <-  seq(nrow(ES_metabolites),1,-1)
  ES_metabolites$score <- (ES_metabolites$score - 1)/(nrow(ES_metabolites) - 1) * 2 + 1
  ES_met_function$score <-  seq(nrow(ES_met_function),1,-1)
  ES_met_function$score <- (ES_met_function$score - 1)/(nrow(ES_met_function) - 1) * 2 + 1
  ES_microbes$score <-  seq(nrow(ES_microbes),1,-1)
  ES_microbes$score <- (ES_microbes$score - 1)/(nrow(ES_microbes) - 1) * 2 + 1
  ES_mic_function$score <-  seq(nrow(ES_mic_function),1,-1)
  ES_mic_function$score <- (ES_mic_function$score - 1)/(nrow(ES_mic_function) - 1) * 2 + 1
}



# III. Calculate abundance score(AS)
{
  metabolites_abun <- met_abun[,DS_metabolites$Node]
  metabolites_function_abun <- met_abun[,DS_met_function$Node]
  AS_metabolites <- apply(metabolites_abun,2,mean)
  AS_metabolites <- data.frame(Node = colnames(metabolites_abun), abundance = AS_metabolites)
  AS_met_function <- apply(metabolites_function_abun,2,mean)
  AS_met_function <- data.frame(Node = colnames(metabolites_function_abun), abundance = AS_met_function)
  
  # Convert the score to[1,3]
  AS_metabolites <- AS_metabolites[order(AS_metabolites$abundance,decreasing = T),]
  AS_metabolites$score <-  seq(nrow(AS_metabolites),1,-1)
  AS_metabolites$score <- (AS_metabolites$score - 1)/(nrow(AS_metabolites) - 1) * 2 + 1
  AS_metabolites$Type <- "metabolite"
  
  # Convert the score to[1,3]
  AS_met_function <- AS_met_function[order(AS_met_function$abundance,decreasing = T),]
  AS_met_function$score <-  seq(nrow(AS_met_function),1,-1)
  AS_met_function$score <- (AS_met_function$score - 1)/(nrow(AS_met_function) - 1) * 2 + 1
  AS_met_function$Type <- "metabolic_function"
  

  microbes_abun <- mic_abun[,DS_microbes$Node]
  microbes_function_abun <- mic_abun[,DS_mic_function$Node]
  AS_microbes <- apply(microbes_abun,2,mean)
  AS_microbes <- data.frame(Node = colnames(microbes_abun), abundance = AS_microbes)
  AS_mic_function <- apply(microbes_function_abun,2,mean)
  AS_mic_function <- data.frame(Node = colnames(microbes_function_abun), abundance = AS_mic_function)
  
  # Convert the score to[1,3]
  AS_microbes <- AS_microbes[order(AS_microbes$abundance,decreasing = T),]
  AS_microbes$score <-  seq(nrow(AS_microbes),1,-1)
  AS_microbes$score <- (AS_microbes$score - 1)/(nrow(AS_microbes) - 1) * 2 + 1
  AS_microbes$Type <- "microbe"
  
  # Convert the score to[1,3]
  AS_mic_function <- AS_mic_function[order(AS_mic_function$abundance,decreasing = T),]
  AS_mic_function$score <-  seq(nrow(AS_mic_function),1,-1)
  AS_mic_function$score <- (AS_mic_function$score - 1)/(nrow(AS_mic_function) - 1) * 2 + 1
  AS_mic_function$group <- "microbial_function"
  
}


# IV. Calculate difference score(DFS)
{
  met_diff$DFS <- sqrt((log10(met_diff$p.adjust))^2 + (log2(met_diff$FC))^2)
  met_diff <- data.frame(Node = met_diff$Node, DFS = met_diff$DFS)
  
  DFS_metabolites <- merge(AS_metabolites,met_diff,by="Node")
  DFS_metabolites$abundance <- DFS_metabolites$DFS
  DFS_metabolites <- DFS_metabolites[,-length(DFS_metabolites)]
  colnames(DFS_metabolites)[2] <- "DFS"
  DFS_metabolites <- DFS_metabolites[order(DFS_metabolites$DFS,decreasing = T),]
  DFS_metabolites$score <-  seq(nrow(DFS_metabolites),1,-1)
  DFS_metabolites$score <- (DFS_metabolites$score - 1)/(nrow(DFS_metabolites) - 1) * 2 + 1
  DFS_metabolites$Type <- "metabolite"
  
  DFS_met_function <- merge(AS_met_function,met_diff,by="Node")
  DFS_met_function$abundance <- DFS_met_function$DFS
  DFS_met_function <- DFS_met_function[,-length(DFS_met_function)]
  colnames(DFS_met_function)[2] <- "DFS"
  DFS_met_function <- DFS_met_function[order(DFS_met_function$DFS,decreasing = T),]
  DFS_met_function$score <-  seq(nrow(DFS_met_function),1,-1)
  DFS_met_function$score <- (DFS_met_function$score - 1)/(nrow(DFS_met_function) - 1) * 2 + 1
  DFS_met_function$Type <- "metabolic_function"
  
  
  mic_diff$DFS <- mic_diff$FC
  mic_diff$DFS <- sqrt((log10(mic_diff$p.value))^2 + (log2(mic_diff$FC))^2)
  mic_diff <- data.frame(Node = mic_diff$Node, DFS = mic_diff$DFS)
  
  DFS_microbes <- merge(AS_microbes,mic_diff,by="Node")
  DFS_microbes$abundance <- DFS_microbes$DFS
  DFS_microbes <- DFS_microbes[,-length(DFS_microbes)]
  colnames(DFS_microbes)[2] <- "DFS"
  DFS_microbes <- DFS_microbes[order(DFS_microbes$DFS,decreasing = T),]
  
  # If FC is inf, then DFS is INF, in which case the score is directly assigned to 3
  flag_inf <- is.finite(DFS_microbes$DFS)
  flag_inf2 <- is.infinite(DFS_microbes$DFS)
  DFS_microbes[flag_inf,]$score <- seq(nrow(DFS_microbes[flag_inf,]),1,-1)
  DFS_microbes[flag_inf,]$score <- (DFS_microbes[flag_inf,]$score - 1)/(nrow(DFS_microbes[flag_inf,]) - 1) * 2 + 1
  DFS_microbes[flag_inf2,]$score <- 3
  DFS_microbes$Type <- "microbe"
  
  
  DFS_mic_function <- merge(AS_mic_function,mic_diff,by="Node")
  DFS_mic_function$abundance <- DFS_mic_function$DFS
  DFS_mic_function <- DFS_mic_function[,-length(DFS_mic_function)]
  colnames(DFS_mic_function)[2] <- "DFS"
  DFS_mic_function <- DFS_mic_function[order(DFS_mic_function$DFS,decreasing = T),]
  
  # If FC is inf, then DFS is INF, in which case the score is directly assigned to 3
  flag_inf <- is.finite(DFS_mic_function$DFS)
  flag_inf2 <- is.infinite(DFS_mic_function$DFS)
  DFS_mic_function[flag_inf,]$score <- seq(nrow(DFS_mic_function[flag_inf,]),1,-1)
  DFS_mic_function[flag_inf,]$score <- (DFS_mic_function[flag_inf,]$score - 1)/(nrow(DFS_mic_function[flag_inf,]) - 1) * 2 + 1
  DFS_mic_function[flag_inf2,]$score <- 3
  DFS_mic_function$Type <- "microbial_function"
}


# IV. Calculate node score(NS)
# The default formula: NS = 0.4*DFS + 0.3*DS + 0.2*ES + 0.1*AS
{
  DS_metabolites <- DS_metabolites[order(DS_metabolites$Node,decreasing = T),]
  DS_met_function <- DS_met_function[order(DS_met_function$Node,decreasing = T),]
  DS_microbes <- DS_microbes[order(DS_microbes$Node,decreasing = T),]
  DS_mic_function <- DS_mic_function[order(DS_mic_function$Node,decreasing = T),]
  
  
  ES_metabolites <- ES_metabolites[order(ES_metabolites$Node,decreasing = T),]
  ES_met_function <- ES_met_function[order(ES_met_function$Node,decreasing = T),]
  ES_microbes <- ES_microbes[order(ES_microbes$Node,decreasing = T),]
  ES_mic_function <- ES_mic_function[order(ES_mic_function$Node,decreasing = T),]
  
  
  AS_metabolites <- AS_metabolites[order(AS_metabolites$Node,decreasing = T),]
  AS_met_function <- AS_met_function[order(AS_met_function$Node,decreasing = T),]
  AS_microbes <- AS_microbes[order(AS_microbes$Node,decreasing = T),]
  AS_mic_function <- AS_mic_function[order(AS_mic_function$Node,decreasing = T),]
  
  
  DFS_metabolites <- DFS_metabolites[order(DFS_metabolites$Node,decreasing = T),]
  DFS_met_function <- DFS_met_function[order(DFS_met_function$Node,decreasing = T),]
  DFS_microbes <- DFS_microbes[order(DFS_microbes$Node,decreasing = T),]
  DFS_mic_function <- DFS_mic_function[order(DFS_mic_function$Node,decreasing = T),]
  
  
  NS_metabolites <- data.frame(Node = DS_metabolites$Node,ns = (0.4*DFS_metabolites$score + 0.3*DS_metabolites$score + 0.2*ES_metabolites$score + 0.1*AS_metabolites$score),Type = DS_metabolites$Type) 
  NS_met_function <- data.frame(Node = DS_met_function$Node,ns = (0.4*DFS_met_function$score + 0.3*DS_met_function$score + 0.2*ES_met_function$score + 0.1*AS_met_function$score),Type = DS_met_function$Type) 
  NS_microbes <- data.frame(Node = DS_microbes$Node,ns = (0.4*DFS_microbes$score + 0.3*DS_microbes$score + 0.2*ES_microbes$score + 0.1*AS_microbes$score),Type = DS_microbes$Type) 
  NS_mic_function <- data.frame(Node = DS_mic_function$Node,ns = (0.4*DFS_mic_function$score + 0.3*DS_mic_function$score + 0.2*ES_mic_function$score + 0.1*AS_mic_function$score),Type = DS_mic_function$Type) 
  
  
  NS_metabolites <- NS_metabolites[order(NS_metabolites$ns,decreasing = T),]
  NS_met_function <- NS_met_function[order(NS_met_function$ns,decreasing = T),]
  NS_microbes <- NS_microbes[order(NS_microbes$ns,decreasing = T),]
  NS_mic_function <- NS_mic_function[order(NS_mic_function$ns,decreasing = T),]
}


write.csv(NS_metabolites,"NS_metabolites.csv",row.names = F)
write.csv(NS_met_function,"NS_met_function.csv",row.names = F)
write.csv(NS_microbes,"NS_microbes.csv",row.names = F)
write.csv(NS_mic_function,"NS_mic_function.csv",row.names = F)