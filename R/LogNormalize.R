LogNormalize <- function(DT, size.factor = 10000){
  DT[,count.normalized:=sum(Count), Cell][]
  DT[,count.normalized:=log1p(Count/count.normalized*size.factor)][]
}

