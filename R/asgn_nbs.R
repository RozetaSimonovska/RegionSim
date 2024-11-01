#' @name asgn_nbs
#' @title Assign neigbors
#'
#' @description Assign k neigbors for regions without neighbors to a neighbors list, based on centroid distance
#'
#' @param neigbs neighbors list
#' @param cen_dist centroid distance matrix
#' @param k number of nearest neigbors to be assigned
#'
#' @return
#' a list of neighbors
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{\link{AggReg}}
#'
#' @examples
#' library("sf")
#' ger <- st_read(system.file(dsn = "shape/gemeinde.shp", package = "RegionSim"))
#' cen_sf <- st_centroid(ger$geometry, of_largest_polygon = TRUE)
#' cen_d <- st_distance(cen_sf)
#' library("units")
#' mun_cen_dist <- drop_units(cen_d)/1000
#' library("spdep")
#' neigbs0 <- poly2nb(ger)
#' neigbs <- asgn_nbs(neigbs = neigbs0, cen_dist = mun_cen_dist, k = 2)
#'
#' @export

asgn_nbs<-function(neigbs, cen_dist, k){
  m<-length(neigbs)
  if(length(neigbs)!=nrow(cen_dist)){
    stop("Length of neibgbours list and distance matrix do not match!")}
  ##### finding regions with no neighbors
  a1<-vector(); for(i in 1:m){a1[i]<-length(neigbs[[i]])}
  a2<-which(a1==1)
  a3<-unlist(neigbs[a2])
  a4<-vector();for(i in which(a3==0)){a4<-c(a4,a2[i])}

  if(length(a4)!=0){
    for(i in a4){
      vec<-c()
      for(j in 1:k){
        vec<-c(vec,which(cen_dist[i,]==cen_dist[i,which(cen_dist[i,]!=0)][
          order(cen_dist[i,which(cen_dist[i,]!=0)])][j]))
      }
      vec2<-vec[1:k]
      vec3<-vec2[order(vec2)]
      if(all(unlist(neigbs[[i]])==0)){
        neigbs[[i]]<-as.vector(vec3)
      }else{ neigbs[[i]]<-as.vector((unique(unlist(c(unlist(neigbs[[i]]),vec3))))) }

     for(j in 1:k){
       if(all(unlist(neigbs[[vec2[j]]]))==0){
         neigbs[[vec2[j]]]<-as.vector(i)
       }else{  neigbs[[vec2[[j]]]]<-as.vector(unique(unlist(c(
         unlist(neigbs[[vec2[j]]]),i)))) }
     }
    }
  }
  return(neigbs)
}
