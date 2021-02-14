#' @name AggReg
#' @title Aggregate regions
#'
#' @description Given \emph{m} smaller regional units aggreagte them into \emph{n} regions based on population restrictions
#'
#' @param sf_pol spatial polygons object
#' @param n the final number of regions for the aggregation
#' @param nseed random seed number
#' @param popvec population vector containig information for the starting m regional units
#' @param areavec area vector containig information for the starting m regional units
#' @param minpop minimum population of the newly created regions. If missing, then set (total population/n)*0.75.
#' @param maxpop maximum population of the newly created regions. If missing, then set (total population/n)*1.25.
#' @param neigbs neighbours list
#' @param cen_dist centroid distances
#' @param n_reg vector of indeces for the regions which will be used to choose the starting n points, default NULL, which means all regional unis from sf_pol will be used
#' @param outpv Variables included in output. Default = FALSE, only provides a spatial polygon of the aggregated regions. If TRUE, additional variables are included (the population and area of the new regions)
#'
#' @return
#' \describe{\emph{newreg}} spatial polygon of the new regions (and list of neighbours \emph{nlist})
#'
#' @author Rozeta Simonovska
#'
#' @import maptools
#' @import rgeos
#'
#' @references
#' Simonovska R., & Tafenau E. (2021). Varying size and shape of spatial units: The MAUP in the case of Germany
#'
#' @seealso \code{\link{asgn_nbs}}
#'
#' @examples
#' library("RegionSim")
#' data(GerMunData, package="RegionSim")
#' library("rgdal")
#' ger<-readOGR(system.file(dsn="shape",package="RegionSim"),layer="gemeinde")
#' ger_lonlat<-spTransform(ger,"+init=epsg:4326")
#' library("rgeos")
#' cen_sf<-gCentroid(ger_lonlat,byid = TRUE)  ### get centroids of the municipalities
#' library("geosphere")
#' sf_mun_cen_dist<-distm(cen_sf@coords)
#' sf_mun_cen_dist<-sf_mun_cen_dist/1000
#' library("spdep")
#' neigbs1<-poly2nb(ger)
#' neigbs<-asgn_nbs(neigbs1,sf_mun_cen_dist,2)
#' ###Connecting the islands
#' ###
#' neigbs[[302]]<-list(unlist(c(unlist(neigbs[[302]]),306)))
#' neigbs[[306]]<-list(unlist(c(unlist(neigbs[[306]]),302)))
#' ###
#' neigbs[[1851]]<-list(unlist(c(unlist(neigbs[[1851]]),1853)))
#' neigbs[[1853]]<-list(unlist(c(unlist(neigbs[[1853]]),1851)))
#' ##
#' neigbs[[9306]]<-list(unlist(c(unlist(neigbs[[9306]]),9395)))
#' neigbs[[9395]]<-list(unlist(c(unlist(neigbs[[9395]]),9306)))
#
#'
#' simn3<-Aggreg(ger,n = 401, popvec = GerMunData$sf_POPULATION, areavec = GerMunData$sf_AREA,
#' minpop = 34400, maxpop = 800000, neigbs = neigbs, cen_dist = sf_mun_cen_dist,
#' n_reg = GerMunData$st_n_reg)
#'
#' @export

Aggreg<-function(sf_pol, n, nseed = 12345, popvec, areavec, minpop=NULL, maxpop=NULL, neigbs, cen_dist, n_reg = NULL, outpv = FALSE ){

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
    large_pop_r<-which(popvec>maxpop)  ##The maximum population
    if(length(large_pop_r)!=0){
      for(im in 1:length(large_pop_r)){
        re[[im]]<-large_pop_r[im]
        n_reg<-n_reg[!n_reg %in% c(re[[im]])]
        n_reg<-n_reg[!n_reg %in% neigbs[[re[[im]]]]]
      }
    }

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
              n_reg2<-n_reg2[!n_reg2 %in% neigbs[[re[[i]]]]]
              if(nseed<start_ns+4){   ###as we increase the seed if it is increased by 3 or less use 80% of average distance else 50%
                  if(n<30 & m/n<3){prcnt<-0.701}else{ prcnt<-0.801 }
              }else{ prcnt<-0.501 }
              n_reg2<-n_reg2[!n_reg2 %in% which(cen_dist[re[[i]],]<sqrt(tot_ar/n)*(prcnt))]
          }
      }
      nseed<-nseed+1
      if(length(re)==n){break}
    }

    nseed<-nseed-1

    ###change order of one of the very large region with another one
    if(length(large_pop_r)!=0){
      r_p1<-sample(1:length(large_pop_r),1)
      r_p2<-sample((length(large_pop_r)+1):n,1)
      pom0<-re[[r_p1]]
      re[[r_p1]]<-re[[r_p2]]
      re[[r_p2]]<-pom0
    }


    frstm<-unlist(re)

    neigbs2<-neigbs    ###availibale neigbours which are not already added to a region
    for (k in 1:m){  neigbs2[[k]]<-list(unlist(neigbs[[k]])[!unlist(neigbs[[k]]) %in% frstm])}
    ###################################################
    ar_nreg<-rep(NA,n)  ###area of the new regions
    for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
    pop_nreg<-rep(NA,n)  ### population of the new region
    for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
    num_neb<-rep(list(NA),n); no_neb_dum<-rep(0,n); all_neb<-vector()  ###no neigbours dummy (indicating that there are no more neigbours to be added to the region)
    for(i in 1:n) {
      for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
      all_neb[i]<-sum(num_neb[[i]])
      if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}


    set.seed(nseed)
    ### adding regions based on population and (having) neigbours
    if(m>5000 & (m/n<50)){
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

    if(m>5000 & (m/n<50)){ ##just samll regions
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
    } ##end if just small regions

    if(no_neb_dum[which(pop_nreg==min(pop_nreg))]==0){  ##If the region with min population still has available neighbours (for the case of small regions)
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
                if(m>5000 & (m/n<50)){
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
    ##############################################
    no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n); all_neb<-vector()
    for(i in 1:n) {
      for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
      all_neb[i]<-sum(num_neb[[i]])
      if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

    for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
    for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
    ################################

    nseed<-nseed+1

    if(pop_nreg[which(pop_nreg==min(pop_nreg))]>min_n_pop || no_neb_dum[which(pop_nreg==min(pop_nreg))]==0){break}

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

  if(m>5000 & (m/n<50) & avrpop/2>minpop){ ##just samll regions
    repeat{
      ###order based on number of neigbours and population, starting with the region with the smallest population number
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
    ###order based on number of neigbours and population, starting with the region with the smallest number
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
              if(m>5000 & (m/n<50)){
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
  } #####END repeat tot_pop/n

  ##############################################
  no_neb_dum<-rep(0,n);num_neb<-rep(list(NA),n)
  for(i in 1:n) {
    for(kk in 1:length(re[[i]])){num_neb[[i]][kk]<-length(unlist(neigbs2[[re[[i]][kk]]]))}
    all_neb[i]<-sum(num_neb[[i]])
    if(all(unlist(num_neb[[i]])==0)){no_neb_dum[i]<-1}}

  for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[re[[i]]]))}
  for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[re[[i]]]))}
  ################################


  ### adding regions based on population, (having) neigbours and area
  set.seed(nseed)
  ##small and medium regions
  if(m>250){
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

    set.seed(nseed)
    ##only small regions
    if(m>5000 & (m/n<50)){
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
              if(m>5000 & (m/n<50)){
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


  ###### Adding the remaining areas to the newly created regions
  ### Create identifier for regions based on letters
  bukvi<-vector()
  for(bi in 1:(length(letters)+1)){for(bj in 1:(length(LETTERS))) {bukvi[bj+length(letters)*(bi-1)]<-paste0(letters[bi-1],LETTERS[bj])}}
  if(n>length(bukvi)){
    lbuk<-floor(n/length(bukvi))+1
    bukvi2<-rep(NA,length(bukvi)*lbuk)
    for(bi2 in 1:(length(bukvi))){
      for(bj2 in 1:lbuk) {
        bukvi2[bi2+length(bukvi)*(bj2-1)]<-paste0(bukvi[bi2],as.character(bj2-1))
      }
    }
    bukvi<-bukvi2
  }
  bukvi<-bukvi[order(bukvi)]

  for(i in 1:n){   sf_Agreg[which(sf_ID %in% re[[i]])]<-bukvi[i]}

  ### Adding the non-assained regions
  set.seed(nseed)
  if(length(unlist(re))!=m){
    repeat{
      no_reg<-which(!sf_Agreg %in% bukvi[1:n])
      if(length(no_reg)<1)  { break
        } else{
          pvar2<-list(); pvar3<-list(); pvar4<-vector()
          for(kk in 1:length(no_reg)){
            pvar2[[kk]]<-unique(unlist(neigbs[[no_reg[kk]]])[!unlist(neigbs[[no_reg[kk]]]) %in% unlist(neigbs2[[no_reg[kk]]])])
            pvar3[[kk]]<-unique(unlist(pvar2[[kk]][!unlist(pvar2[[kk]]) %in% no_reg]))
            pvar4[kk]<-length(unlist(unlist(pvar3[[kk]])))
          }
          pvarl<-pvar3; pvar5<-vector()
          for(ik in which(pvar4!=0)){
            for(il in 1:length(pvar3[[ik]])){
              if(pop_nreg[which(bukvi==sf_Agreg[unlist(pvar3[[ik]][il])])]>=avrpop){
                pvarl[[ik]]<-pvarl[[ik]][!pvar3[[ik]] %in% pvar3[[ik]][il]]
                pvarl[[ik]]<-pvarl[[ik]][!is.na(pvarl[[ik]])]
              }}
          }
          for(kk in 1:length(no_reg)){pvar5[kk]<-length(unlist(unlist(pvarl[[kk]])))}

          if(length(which(pvar5!=0))>0){
            pvar6<-as.numeric(sample(as.character(which(pvar5!=0)),1))
            pom4<-as.numeric(sample(as.character(pvarl[[pvar6]]),1))
            } else{
            pvar6<-as.numeric(sample(as.character(which(pvar4!=0)),1))
            pom5<-vector()
            for(im in 1:length(pvar3[[pvar6]])){  pom5[im]<-pop_nreg[which(bukvi==sf_Agreg[pvar3[[pvar6]][im]])]}
            pom4<-pvar3[[pvar6]][which(pom5==min(pom5))][1]
          }

          sf_Agreg[no_reg[pvar6]]<-sf_Agreg[pom4]
          for (k in 1:m){ neigbs2[[k]]<-list(unlist(neigbs2[[k]])[!unlist(neigbs2[[k]]) %in% no_reg[pvar6]])}
      }
    }###END Repeat
  }


  ar_nreg<-rep(NA,n)
  for (i in 1:n){ar_nreg[i]<-as.numeric(sum(areavec[which(sf_Agreg==bukvi[i])]))}
  pop_nreg<-rep(NA,n)
  for (i in 1:n){pop_nreg[i]<-as.numeric(sum(popvec[which(sf_Agreg==bukvi[i])]))}
  ###################################################################################
  ###################################################
  ###############################
  new_reg_id<-vector()
  for(ij in 1:m) {new_reg_id[ij]<-which(sf_Agreg[ij]==bukvi)}

  poly <- rgeos::gBuffer(sf_pol, byid=TRUE, width=0)
  newreg<-maptools::unionSpatialPolygons(poly,new_reg_id)


  ######################################################
  if(outpv){
    varbl<-list(ar_nreg,pop_nreg,new_reg_id,nseed); names(varbl)<-c("sf_ar","sf_pop","new_reg_id","nseed")
    all_rez<-list(newreg,varbl)
    names(all_rez)<-c("newreg","varbl")

    return(all_rez)
  } else{ return(newreg)}


}
