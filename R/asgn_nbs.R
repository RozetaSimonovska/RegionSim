#' @name asgn_nbs
#' @title Assign neigbours
#'
#' @description Assign k neigbours for regions without neighbours to a neighbours list, based on centroid distance
#'
#' @param neigbs neighbours list
#' @param cen_dist centroid distance matrix
#' @param k number of nearest neigbours to be assigned
#'
#' @return
#' a list of neighbours
#'
#' @author Rozeta Simonovska
#'
#' @references
#' Simonovska R., & Tafenau E. (2021). Varying size and shape of spatial units: The MAUP in the case of Germany
#'
#' @seealso \code{\link{AggReg}}
#'
#' @export

asgn_nbs<-function(neigbs, cen_dist, k){
  m<-length(neigbs)
  if(length(neigbs)!=nrow(cen_dist)){stop("Length of neibgbours list and distance matrix do not match!")}
  ##### finding regions with no neigbours
  a1<-vector(); for(i in 1:m){a1[i]<-length(neigbs[[i]])}
  a2<-unlist(neigbs[which(a1==1)]); a3<-which(a1==1)
  a4<-vector();for(i in which(a2==0)){a4<-c(a4,a3[i])}

  if(length(a4)!=0){
    for(i in a4){
      vec<-c()
      for(j in 1:k){
        vec<-c(vec,which(cen_dist[i,]==cen_dist[i,which(cen_dist[i,]!=0)][order(cen_dist[i,which(cen_dist[i,]!=0)])][j]))
      }
      vec2<-vec[1:k]
      vec3<-vec2[order(vec2)]
      if(all(unlist(neigbs[[i]])==0)){
        neigbs[[i]]<-list(c(vec3))
      }else{ neigbs[[i]]<-list(unique(unlist(c(unlist(neigbs[[i]]),vec3)))) }

     for(j in 1:k){
       if(all(unlist(neigbs[[vec2[j]]]))==0){
         neigbs[[vec2[j]]]<-list(c(i))
       }else{  neigbs[[vec2[[j]]]]<-list(unique(unlist(c(unlist(neigbs[[vec2[j]]]),i)))) }
     }
    }
  }
  return(neigbs)
}
