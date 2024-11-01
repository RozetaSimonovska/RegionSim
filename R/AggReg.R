#' @name AggReg
#' @title Aggregate regions
#'
#' @description Given \emph{m} smaller regional units aggregate them into \emph{n} regions based on population restrictions
#'
#' @param sf_pol spatial polygons object
#' @param n the final number of regions for the aggregation
#' @param nseed random seed number
#' @param popvec population vector containing information for the starting m regional units
#' @param areavec area vector containing information for the starting m regional units
#' @param minpop minimum population of the newly created regions. If missing, then set (total population/n)*0.75.
#' @param maxpop maximum population of the newly created regions. If missing, then set (total population/n)*1.25.
#' @param neigbs neighbors list
#' @param cen_dist centroid distances
#' @param n_reg vector of indexes for the regions which will be used to choose the starting n points, default NULL, which means all regional unis from sf_pol will be used
#' @param outpv variables included in output. Default = FALSE, only provides a spatial polygon of the aggregated regions. If TRUE, additional variables are included (the population and area of the new regions)
#' @param nopoprist Logical. If there should be no population restrictions
#' @param onebyone Logical. I remaining regional units at the end should be added one by one
#' @param repuntilmin Logical. If for minimum population condition has to be fulfilled for all newly created regions
#'
#' @return
#' \describe{\emph{newreg}} spatial polygon of the new regions (and list of neighbors \emph{nlist})
#'
#' @author Rozeta Simonovska
#'
#' @importFrom spdep poly2nb
#' @import sf
#' @importFrom units drop_units
#' @import dplyr
#'
#' @references
#' Simonovska R., & Tafenau E. (2024). Varying size and shape of spatial units: Analysing the MAUP through agglomeration economies in the case of Germany. REGION, 11(2), 63–97.
#'
#' @seealso \code{\link{asgn_nbs}}
#'
#' @examples
#' library("RegionSim")
#' library("sf")
#' gerNUTS3<-st_read(system.file(dsn = "shape/kreise.shp", package = "RegionSim"))
#'
#' simN2<-Aggreg(gerNUTS3,
#'               n = 38,
#'               popvec = gerNUTS3$popNUTS3,
#'               areavec = gerNUTS3$areaNUTS3,
#'               minpop = 678750,
#'               maxpop = 7000000,
#'               nopoprist = TRUE,
#'               onebyone = FALSE)
#' plot(simN2)
#'
#' @export

Aggreg<-function(sf_pol=NULL, n, nseed = 12345,
                 popvec, areavec, minpop=NULL, maxpop=NULL,
                 neigbs=NULL, cen_dist=NULL, n_reg = NULL, outpv = FALSE,
                 nopoprist = FALSE,
                 onebyone = TRUE,
                 repuntilmin = FALSE
                 ){

  if(is.null(sf_pol) & is.null(neigbs) & is.null(cen_dist)){
    stop("Missing sf_pol, neigbs and cen_dist")
  }else {
    if(is.null(neigbs) ){
      neigbs <- spdep::poly2nb(sf_pol)
    }
    if(is.null(cen_dist)){
      cen_sf <- sf::st_centroid(st_geometry(sf_pol), of_largest_polygon = TRUE)
      cen_d <- sf::st_distance(cen_sf)
      cen_dist <- units::drop_units(cen_d)/1000
    }
  }
  m <- length(neigbs)
  if(length(neigbs)!=nrow(cen_dist)){stop("Length of neibgs and sf_dist do not match!")}

  tot_ar <- sum(areavec); tot_pop <- sum(popvec)
  avrpop <- tot_pop/n; avrarea <- tot_ar/n

  if(is.null(minpop)){ minpop <- avrpop*0.75}
  if(is.null(maxpop)){ maxpop <- avrpop*1.25}

  if(is.null(n_reg)){ n_reg<-seq(1,m,1) }

  sf_ID<-seq(1,m,1); sf_Agreg<-seq(1,m,1)

  start_ns<-nseed ###starting random seed

  repeat{ ### repeat unit minimal population in each new region is over minpop
    set.seed(nseed)  ##set seed
    re<-list()

    st_seed2<-nseed
    ##############
    #add regions with over average population
    if(!nopoprist){
      large_pop_r<-which(popvec>maxpop)  ##The maximum population
      if(length(large_pop_r)!=0){
        for(im in 1:length(large_pop_r)){
          re[[im]]<-large_pop_r[im]
          n_reg<-n_reg[!n_reg %in% c(re[[im]])]
          n_reg<-n_reg[!n_reg %in% neigbs[[re[[im]]]]]
        }
      }
    }else{  large_pop_r<-vector()  }

    nseed<-nseed-1
    ###add remaining starting regions
    repeat{
      nseed<-nseed+1
      set.seed(nseed)
      n_reg2<-n_reg
      for(i in (length(large_pop_r)+1):n){
        if(length(n_reg2)==0){break
        } else {
          re[[i]]<-as.numeric(sample(as.character(n_reg2),1))
          n_reg2<-n_reg2[!n_reg2 %in% c(re[[i]])]
          if(!(n<30 & m/n<3)){n_reg2<-n_reg2[!n_reg2 %in% neigbs[[re[[i]]]]]}
          if(nseed<start_ns+4){   ###as we increase the seed if it is increased by 3 or less use 80% of average distance else 50%
            if(n<30 & m/n<3){prcnt<-0.701}else{ prcnt<-0.801 }
          }else{ prcnt<-0.501 }
          n_reg2<-n_reg2[!n_reg2 %in% which(cen_dist[re[[i]],]<sqrt(tot_ar/n)*(prcnt))]
        }
      }
      if(length(re)==n){break}
    }

    ##nseed<-nseed-1

    ###change order of one of the very large region with another one
    if(!nopoprist){
      if(length(large_pop_r)!=0){
        r_p1<-sample(1:length(large_pop_r),1)
        r_p2<-sample((length(large_pop_r)+1):n,1)
        temp0<-re[[r_p1]]
        re[[r_p1]]<-re[[r_p2]]
        re[[r_p2]]<-temp0
      }
    }


    frstm<-unlist(re)

    neigbs2<-neigbs    ###available neighbors which are not already added to a region
    for (k in 1:m){  neigbs2[[k]]<-list(unlist(neigbs[[k]])[!unlist(neigbs[[k]]) %in% frstm])}
    ###################################################
    ar_nreg<-rep(NA,n)  ###area of the new regions
    for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
    pop_nreg<-rep(NA,n)  ### population of the new region
    for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
    num_neb<-rep(list(NA),n); no_neb_dum<-rep(0,n); all_neb<-vector()  ###no neighbors dummy (indicating that there are no more neighbors to be added to the region)
    for(i in 1:n) {
      for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
      all_neb[i]<-sum(num_neb[[i]])
      if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}


    set.seed(nseed)
    ### adding regions based on population and (having) neighbors
    if((m>5000 & (m/n<50)) || repuntilmin){
      if(nseed<(st_seed2+2)){
        min_n_pop<-minpop
      }else if(nseed<(st_seed2+4)){
        min_n_pop<-minpop*0.8
      }else if(nseed<(st_seed2+6)){
        min_n_pop<-minpop*0.66667
      }else{
        min_n_pop<-minpop*0.5
      }
    }else{ min_n_pop<-minpop }

    if(!nopoprist){
      if((m>5000 & (m/n<50)) || repuntilmin){ ##just small regions
        repeat{
          ###order based on population, starting with the region with the smallest number
          for(i in which(pop_nreg<min_n_pop & no_neb_dum==0)[
            order(pop_nreg[which(pop_nreg<min_n_pop & no_neb_dum==0)],
                  all_neb[which(pop_nreg<min_n_pop & no_neb_dum==0)])]) {
            if(pop_nreg[i]>min_n_pop || no_neb_dum[i]==1){ next
            } else{
              neb_l<-list()
              for(kk in 1:length(re[[i]])){ neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])  }
              neb_v<-unique(unlist(neb_l))
              if(length(neb_v)==0){ next
              } else{
                pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
                re[[i]]<-c(re[[i]],pvar1)
                for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1]) }
                ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
                pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
                no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
                for(i in 1:n) {
                  for(kk in 1:length(re[[i]])){ num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]])) }
                  all_neb[i]<-sum(num_neb[[i]])
                  if(all(unlist(num_neb[[i]])==0)){ no_neb_dum[i]<-1 }
                }
              }
            }
          }
          if(length( which(pop_nreg<min_n_pop & no_neb_dum==0))==0)     {break}
        }##end repeat

        ##############################################
        no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n); all_neb<-vector()
        for(i in 1:n) {
          for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
          all_neb[i]<-sum(num_neb[[i]])
          if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

        for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
        for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
        ################################


        if(any(no_neb_dum[which(pop_nreg==min(pop_nreg))]==0)){  ##If the region with min population still has available neighbors (for the case of small regions)
          repeat{
            ###order based on population, starting with the region with the smallest number
            for(i in which(pop_nreg<minpop & no_neb_dum==0)[
              order(all_neb[which(pop_nreg<minpop & no_neb_dum==0)],
                    pop_nreg[which(pop_nreg<minpop & no_neb_dum==0)])]) {
              if(pop_nreg[i]>minpop || no_neb_dum[i]==1){ next
              } else{
                neb_l<-list()
                for(kk in 1:length(re[[i]])){ neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])  }
                neb_v<-unique(unlist(neb_l))
                if(length(neb_v)==0){ next
                }  else{
                  if((m>5000 & (m/n<50)) || repuntilmin){
                    pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
                  }else { pvar1<-as.numeric(sample(as.character(neb_v),1))  }
                  re[[i]]<-c(re[[i]],pvar1)
                  for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
                  ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
                  pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
                  no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
                  for(i in 1:n) {
                    for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
                    all_neb[i]<-sum(num_neb[[i]])
                    if(all(unlist(num_neb[[i]])==0)){  no_neb_dum[i]<-1  }
                  }
                }
              }
            }
            if(length( which(pop_nreg<minpop & no_neb_dum==0))==0 || (
              pop_nreg[which(pop_nreg==min(pop_nreg))]<min_n_pop & no_neb_dum[which(pop_nreg==min(pop_nreg))]==1)){break}

          }##end repeat
        }###end if
      } ##end if just small regions
      ##############################################
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n); all_neb<-vector()
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
      for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
      ################################
    }
    nseed<-nseed+1

    if(all(pop_nreg[which(pop_nreg==min(pop_nreg))]>min_n_pop) ||
       all(no_neb_dum[which(pop_nreg==min(pop_nreg))]==0) || nopoprist){break}

  } ### END repeat unit minimal population in each new region is over minpop (start =starting regions)
  #######################################################################################################
  nseed<-nseed-1
  set.seed(nseed)

  ##############################################
  no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
  for(i in 1:n) {
    for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
    all_neb[i]<-sum(num_neb[[i]])
    if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

  for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
  for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}

  if(!nopoprist){
    if(((m>5000 & (m/n<50)) || repuntilmin) & avrpop/2>minpop){ ##just small regions
      repeat{
        ###order based on number of neighbors and population, starting with the region with the smallest population number
        for(i in which(pop_nreg<avrpop/2 & no_neb_dum==0)[
          order(all_neb[which(pop_nreg<(avrpop)/2 & no_neb_dum==0)],
                pop_nreg[which(pop_nreg<(avrpop)/2 & no_neb_dum==0)])]) {
          if(pop_nreg[i]>avrpop/2 || no_neb_dum[i]==1){ next
          } else{
            neb_l<-list()
            for(kk in 1:length(re[[i]])){ neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])  }
            neb_v<-unique(unlist(neb_l))
            if(length(neb_v)==0){ next
            } else{
              pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
              re[[i]]<-c(re[[i]],pvar1)
              for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
              ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
              pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
              no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
              for(i in 1:n) {
                for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
                all_neb[i]<-sum(num_neb[[i]])
                if(all(unlist(num_neb[[i]])==0)){  no_neb_dum[i]<-1  }
              }
            }
          }
        }
        ####
        if(length(which(pop_nreg<avrpop/2 & no_neb_dum==0))==0) {break}
      } ####END repeat avrpop/2

      ##############################################
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
      for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
      ################################
    }###end if ## just small regions

    repeat{
      ###order based on number of neighbors and population, starting with the region with the smallest number
      for(i in which(pop_nreg<avrpop & no_neb_dum==0)[
        order(all_neb[which(pop_nreg<avrpop & no_neb_dum==0)],
              pop_nreg[which(pop_nreg<avrpop & no_neb_dum==0)])]) {
        if(pop_nreg[i]>avrpop || no_neb_dum[i]==1){ next
        } else{
          neb_l<-list()
          for(kk in 1:length(re[[i]])){ neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])  }
          neb_v<-unique(unlist(neb_l))
          if(length(neb_v)==0){ next
          } else{
            if((m>5000 & (m/n<50)) || repuntilmin){
              pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
            }else { pvar1<-as.numeric(sample(as.character(neb_v),1))}
            re[[i]]<-c(re[[i]],pvar1)
            for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
            ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
            pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
            no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
            for(i in 1:n) {
              for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
              all_neb[i]<-sum(num_neb[[i]])
              if(all(unlist(num_neb[[i]])==0)){ no_neb_dum[i]<-1 }
            }
          }
        }
      }
      ####
      if(length(which(pop_nreg<avrpop & no_neb_dum==0))==0)     {break}
    } #####END repeat avrpop

    ##############################################
    no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
    for(i in 1:n) {
      for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
      all_neb[i]<-sum(num_neb[[i]])
      if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

    for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
    for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
    # ################################


    ### adding regions based on population, (having) neighbors and area
    set.seed(nseed)
    ##small and medium regions
    if(m>250){
      if(length(which(ar_nreg<avrarea & pop_nreg<avrpop & no_neb_dum==0)[
        order(pop_nreg[which(ar_nreg<avrarea & pop_nreg<avrpop & no_neb_dum==0)])])>0){
        repeat{
          ###order based on population, starting with the region with the smallest number
          for(i in which(ar_nreg<avrarea & pop_nreg<avrpop & no_neb_dum==0)[
            order(pop_nreg[which(ar_nreg<avrarea & pop_nreg<avrpop & no_neb_dum==0)])]){
            if(ar_nreg[i]>avrarea || pop_nreg[i]>avrpop || no_neb_dum[i]==1){next
            } else{
              neb_l<-list()
              for(kk in 1:length(re[[i]])){  neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])   }
              neb_v<-unique(unlist(neb_l))
              if(length(neb_v)==0){ next
              } else{
                pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
                re[[i]]<-c(re[[i]],pvar1)
                for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
                ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
                pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
                no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
                for(i in 1:n) {
                  for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
                  all_neb[i]<-sum(num_neb[[i]])
                  if(all(unlist(num_neb[[i]])==0)){ no_neb_dum[i]<-1 }
                }
              }
            }
          }
          no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
          for(i in 1:n) {
            for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
            all_neb[i]<-sum(num_neb[[i]])
            if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

          if(length(which(ar_nreg<avrarea & pop_nreg<avrpop*2 & no_neb_dum==0))==0){break}
        } ###END repeat ar_nreg<avrarea & pop_nreg<avrpop
      }##end if

      ########################
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
      for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
      ###############################
    }####END if m>250 ##
  }

  ###no rist
  if(nopoprist){
    repeat{
      ###order based on population, starting with the region with the smallest number
      for(i in which(ar_nreg<avrarea*0.7 & no_neb_dum==0)){
        if(ar_nreg[i]>avrarea*0.7 || no_neb_dum[i]==1){next
        } else{
          neb_l<-list()
          for(kk in 1:length(re[[i]])){  neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])   }
          neb_v<-unique(unlist(neb_l))
          if(length(neb_v)==0){ next
          } else{
            if(m/n>250){
              if(length(neb_v)>=100){sizen <- 100
              }else{sizen <- length(neb_v)}
            }else if(m/n>100){
              if(length(neb_v)>=50){sizen <- 50
              }else{sizen <- length(neb_v)}
            }else if(m/n>15){
              if(length(neb_v)>=10){sizen <- 10
              }else{sizen <- length(neb_v)}
            }else{sizen <- length(neb_v)}
            pvar1<-as.numeric(sample(as.character(neb_v),sizen))
            re[[i]]<-c(re[[i]],pvar1)
            for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
            ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
            pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
            no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
            for(i in 1:n) {
              for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
              all_neb[i]<-sum(num_neb[[i]])
              if(all(unlist(num_neb[[i]])==0)){ no_neb_dum[i]<-1 }
            }
          }
        }
      }
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      if(length(which(ar_nreg<avrarea*0.7 & no_neb_dum==0))==0){break}
    } ###END repeat ar_nreg<avrarea*0.7
    no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
    for(i in 1:n) {
      for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
      all_neb[i]<-sum(num_neb[[i]])
      if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

    for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
    for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}

    repeat{
      ###order based on population, starting with the region with the smallest number
      for(i in which(ar_nreg<avrarea*0.9 & no_neb_dum==0)){
        if(ar_nreg[i]>avrarea*0.9 || no_neb_dum[i]==1){next
        } else{
          neb_l<-list()
          for(kk in 1:length(re[[i]])){  neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])   }
          neb_v<-unique(unlist(neb_l))
          if(length(neb_v)==0){ next
          } else{
            if(length(neb_v)>=5){sizen <- 5
            }else{sizen <- length(neb_v)}
            pvar1<-as.numeric(sample(as.character(neb_v),sizen))
            re[[i]]<-c(re[[i]],pvar1)
            for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
            ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
            pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
            no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
            for(i in 1:n) {
              for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
              all_neb[i]<-sum(num_neb[[i]])
              if(all(unlist(num_neb[[i]])==0)){ no_neb_dum[i]<-1 }
            }
          }
        }
      }
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      if(length(which(ar_nreg<avrarea*0.9 & no_neb_dum==0))==0){break}
    } ###END repeat ar_nreg<avrarea*0.9
  }##no rist
  no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
  for(i in 1:n) {
    for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
    all_neb[i]<-sum(num_neb[[i]])
    if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

  for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
  for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}

  set.seed(nseed)
  if(!nopoprist){
    ##only small regions
    if((m>5000 & (m/n<50)) || repuntilmin){
      repeat{
        ###order based on population, starting with the region with the smallest number
        for(i in which(ar_nreg<avrarea*2 & pop_nreg<avrpop*2 & no_neb_dum==0)[
          order(pop_nreg[which(ar_nreg<avrarea*2 & pop_nreg<avrpop*2 & no_neb_dum==0)])]){
          neb_l<-list()
          for(kk in 1:length(re[[i]])){  neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])   }
          neb_v<-unique(unlist(neb_l))
          if(length(neb_v)==0){ next
          } else{
            pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
            re[[i]]<-c(re[[i]],pvar1)
            for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
            ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
            pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
            no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
            for(i in 1:n) {
              for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
              all_neb[i]<-sum(num_neb[[i]])
              if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}
            }
          }
        }
        no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
        for(i in 1:n) {
          for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
          all_neb[i]<-sum(num_neb[[i]])
          if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

        if(length(which(ar_nreg<avrarea*2 & pop_nreg<avrpop*2 & no_neb_dum==0))==0){break}
      } ####END repeat avrarea*2 & pop_nreg<avrpop*2 #

      ########################
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
      for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}

      ###############################
      set.seed(nseed)
      repeat{
        ###order based on population, starting with the region with the smallest number
        for(i in which(ar_nreg<avrarea*2 & pop_nreg<maxpop & no_neb_dum==0)[
          order(pop_nreg[which(ar_nreg<avrarea*2 & pop_nreg<maxpop & no_neb_dum==0)])]){
          neb_l<-list()
          for(kk in 1:length(re[[i]])){  neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])   }
          neb_v<-unique(unlist(neb_l))
          if(length(neb_v)==0){ next
          } else{
            pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
            re[[i]]<-c(re[[i]],pvar1)
            for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
            ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
            pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
            no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
            for(i in 1:n) {
              for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
              all_neb[i]<-sum(num_neb[[i]])
              if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}
            }
          }
        }
        no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
        for(i in 1:n) {
          for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
          all_neb[i]<-sum(num_neb[[i]])
          if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

        if(length(which(ar_nreg<avrarea*2 & pop_nreg<maxpop & no_neb_dum==0))==0){break}
      }###END repeat  ar_nreg<avrarea*2 & pop_nreg<maxpop #

      ########################
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
      for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
    }

    ############################################################################
    set.seed(nseed)
    if(length(unlist(re))!=m){
      repeat{
        ###order based on population, starting with the region with the smallest number
        for(i in which(pop_nreg<maxpop & no_neb_dum==0)[
          order(all_neb[which(pop_nreg<maxpop & no_neb_dum==0)],
                pop_nreg[which(pop_nreg<maxpop & no_neb_dum==0)])]) {
          if(pop_nreg[i]>maxpop || no_neb_dum[i]==1){next
          } else{
            neb_l<-list()
            for(kk in 1:length(re[[i]])){  neb_l[[kk]]<-unlist(neigbs2[[re[[i]][kk]]])   }
            neb_v<-unique(unlist(neb_l))
            if(length(neb_v)==0){ next
            } else{
              if((m>5000 & (m/n<50)) || repuntilmin){
                pvar1<-neb_v[which(popvec[neb_v]==max(popvec[neb_v]))]
              }else { pvar1<-as.numeric(sample(as.character(neb_v),1))}
              re[[i]]<-c(re[[i]],pvar1)
              for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% pvar1])}
              ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))
              pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))
              no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
              for(i in 1:n) {
                for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
                all_neb[i]<-sum(num_neb[[i]])
                if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}
              }
            }
          }
        }
        no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
        for(i in 1:n) {
          for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
          all_neb[i]<-sum(num_neb[[i]])
          if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

        if(length(which(pop_nreg<maxpop & no_neb_dum==0))==0){break}
      } ###END repeat pop_nreg<maxpop #

      ########################
      no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
      for(i in 1:n) {
        for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
        all_neb[i]<-sum(num_neb[[i]])
        if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

      for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
      for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
    }###END IF
  }


  ###### Adding the remaining areas to the newly created regions
  ### Create identifier for regions based on letters
  letters_vec<-vector()
  for(bi in 1:(length(letters)+1)){for(bj in 1:(length(LETTERS))) {letters_vec[bj+length(letters)*(bi-1)]<-paste0(letters[bi-1],LETTERS[bj])}}
  if(n>length(letters_vec)){
    lbuk<-floor(n/length(letters_vec))+1
    letters_vec2<-rep(NA,length(letters_vec)*lbuk)
    for(bi2 in 1:(length(letters_vec))){
      for(bj2 in 1:lbuk) {
        letters_vec2[bi2+length(letters_vec)*(bj2-1)]<-paste0(letters_vec[bi2],as.character(bj2-1))
      }
    }
    letters_vec<-letters_vec2
  }
  letters_vec<-letters_vec[order(letters_vec)]

  for(i in 1:n){   sf_Agreg[which(sf_ID %in% re[[i]])]<-letters_vec[i]}

  ### Adding the non-assigned regions
  set.seed(nseed)
  if(length(unlist(re))!=m){
    if(onebyone){
      repeat{
        no_reg<-which(!sf_Agreg %in% letters_vec[1:n])
        if(length(no_reg)>=1)  {
          pvar2<-list(); pvar3<-list(); pvar4<-vector()
          for(kk in 1:length(no_reg)){

            pvar2[[kk]]<-unique(unlist(neigbs[[no_reg[kk]]])[!unlist(neigbs[[no_reg[kk]]]) %in% unlist(neigbs2[[no_reg[kk]]])])

            pvar3[[kk]]<-unique(unlist(pvar2[[kk]][!unlist(pvar2[[kk]]) %in% no_reg]))

            pvar4[kk]<-length(unlist(unlist(pvar3[[kk]])))

          }
          pvarl<-pvar3; pvar5<-vector()
          if(!nopoprist){
            for(ik in which(pvar4!=0)){
              for(il in 1:length(pvar3[[ik]])){
                if(pop_nreg[which(letters_vec==sf_Agreg[unlist(pvar3[[ik]][il])])]>=avrpop){
                  pvarl[[ik]]<-pvarl[[ik]][!pvar3[[ik]] %in% pvar3[[ik]][il]]
                  pvarl[[ik]]<-pvarl[[ik]][!is.na(pvarl[[ik]])]
                }}
            }
          }
          for(kk in 1:length(no_reg)){pvar5[kk]<-length(unlist(unlist(pvarl[[kk]])))}

          if(length(which(pvar5!=0))>0){
            pvar6<-as.numeric(sample(as.character(which(pvar5!=0)),1))
            temp4<-as.numeric(sample(as.character(pvarl[[pvar6]]),1))
          } else{
            pvar6<-as.numeric(sample(as.character(which(pvar4!=0)),1))
            temp5<-vector()
            if(!nopoprist){
              for(im in 1:length(pvar3[[pvar6]])){  temp5[im]<-pop_nreg[which(letters_vec==sf_Agreg[pvar3[[pvar6]][im]])]}
              temp4<-pvar3[[pvar6]][which(temp5==min(temp5))][1]
            }else {
              temp4<-pvar3[[pvar6]][1]
            }
          }

          sf_Agreg[no_reg[pvar6]]<-sf_Agreg[temp4]
          for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% no_reg[pvar6]])}
        }
        if(length(no_reg)<1)  {break}
      }###END Repeat
    }else{
      nind<-0
      repeat{
        nind<-nind+1
        no_reg<-which(!sf_Agreg %in% letters_vec[1:n])
        #print(nind); print(length(no_reg))
        if(length(no_reg)>=1){
          nbslv<-sapply(1:m, function(x) length(unlist(neigbs2[[x]])))
          for(is in which(nbslv!=0)){
            if(length(unlist(neigbs2[[is]]))>0 & (sf_Agreg[is] %in% letters_vec[1:n])){
              ir<-0
              repeat{
                ir <- ir+1
                if((sf_Agreg[is] %in% letters_vec[ir]) || ir==n){ break}
              }
              tmp<-unlist(neigbs2[[is]])[which(unlist(neigbs2[[is]]) %in% no_reg)]
              if(length(tmp)>0){
                sf_Agreg[tmp]<-letters_vec[ir]
                for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% tmp])}
              }
            }
          }
        }
        if(length(no_reg)<1)  {break}
      }
    }
  }

  ar_nreg<-rep(NA,n)
  for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[which(sf_Agreg==letters_vec[i])]))}
  pop_nreg<-rep(NA,n)
  for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[which(sf_Agreg==letters_vec[i])]))}
  ###################################################################################
  ###################################################
  ###############################
  new_reg_id<-vector()
  for(ij in 1:m) {new_reg_id[ij]<-which(sf_Agreg[ij]==letters_vec)}

  poly <-  sf::st_buffer(sf_pol, dist=1000)
  poly$id <- new_reg_id
  newreg<- poly %>% dplyr::select(id) %>%
    dplyr::group_by(id)  %>% dplyr::summarise()



  ######################################################
  if(outpv){
    varbl<-list(ar_nreg, pop_nreg, new_reg_id, nseed)
    names(varbl)<-c("sf_ar","sf_pop","new_reg_id","nseed")

    all_rez<-list(newreg,varbl)
    names(all_rez)<-c("newreg","varbl")

    return(all_rez)
  } else{ return(newreg)}


}

