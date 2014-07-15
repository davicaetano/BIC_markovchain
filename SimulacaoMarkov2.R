library(expm)

matequal <- function(x, y)
    is.matrix(x) && is.matrix(y) && all(dim(x) == dim(y)) && all(x ==y)


SIMULA <- function(A,N)
{    
    k = round(log(nrow(A),3)) 
    SAIDA <- sample(0:2,k,replace=TRUE)
    for(i in ((k):(N-1))){
        LINHA = 1
        for (l in 0:(k-1))
        {
            LINHA = LINHA + SAIDA[i-l] * (3 ^ (l))
        }
        SOMA <- 0
        ALE <- runif(1)
        for(j in 1:ncol(A)){        
            SOMA <- SOMA + A[LINHA,j]
            if(ALE <= SOMA){
                SAIDA <- c(SAIDA,j-1)
                ALE <- 2
            }
        }
    }
    return (SAIDA)
}

ESTIMA <- function(VET)
{
    BIC = 1:10 
    for(n in 1:10){ #Estou considerando ordem no mínimo 1 e no máximo 10
        RES <- matrix(rep(0,3^n),nrow=3^n,ncol=3)
        for(i in n:(length(VET)-1))
        {
            LINHA = 1
            for (l in 0:(n-1))
            {
                LINHA = LINHA + VET[i-l] * (3 ^ (l))
            }
            RES[LINHA,VET[i+1]+1] = RES[LINHA,VET[i+1]+1] + 1
        }
        RES2 = RES / RES %*% rep(1,ncol(RES)) %*% t(rep(1,ncol(RES)))
        Like = 0
        for(i in n:(length(VET)-1))
        {
            LINHA = 1
            for (l in 0:(n-1))
            {
                LINHA = LINHA + VET[i-l] * (3 ^ (l))
            }
            Like = Like + log(RES2[LINHA,VET[i+1]+1])
        }
        BIC[n] = -2*Like + (3 ^ (n+1)) * 2 * (3 - 1) * log(n)/2
    }
    RESP = 1:10
    return(RESP[BIC==min(BIC)])
}

RODAR <- function(k,N1,N2=20){
    l = runif(3^k*2)
    A = c(min(l[1],l[2]),max(l[2],l[1])-min(l[1],l[2]),1-max(l[2],l[1]))
    for (i in 1:(3^k-1)){
        A = rbind(A,c(min(l[2*i+1],l[2*i+2]),max(l[2*i+2],l[2*i+1])-min(l[2*i+1],l[2*i+2]),1-max(l[2*i+2],l[2*i+1])))
    }
    B<-SIMULA(A,N1)
    
    PI <- matrix(0,2,nrow=N2)
    colnames(PI) = c("Passos","Ordem estimada")
    for (i in 1:N2)
    {
        PI[i,1] <- i*N1/N2
        PI[i,2] <- ESTIMA(B[1:(i*N1/N2)])
    }
    return(PI)
}

RODAR(3,1000)
RODAR(4,10000)
RODAR(5,10000)
RODAR(6,100000)

ESTIMA2 <- function(VET)
{
    BIC = 1:10 
    for(n in 1:10){
        RES <- matrix(rep(0,3^n),nrow=3^n,ncol=3)
        for(i in n:(length(VET)-1))
        {
            LINHA = 1
            for (l in 0:(n-1))
            {
                LINHA = LINHA + VET[i-l] * (3 ^ (l))
            }
            RES[LINHA,VET[i+1]+1] = RES[LINHA,VET[i+1]+1] + 1
        }
        RES2 = RES / RES %*% rep(1,ncol(RES)) %*% t(rep(1,ncol(RES)))
        Like = 0
        for(i in n:(length(VET)-1))
        {
            LINHA = 1
            for (l in 0:(n-1))
            {
                LINHA = LINHA + VET[i-l] * (3 ^ (l))
            }
            Like = Like + log(RES2[LINHA,VET[i+1]+1])
        }
        BIC[n] = -2*Like + (3 ^ (n+1)) * 2 * (3 - 1) * log(n)/2
    }
    return(BIC)
}

RODAR2 <- function(k,N1,N2=20){
    l = runif(3^k*2)
    A = c(min(l[1],l[2]),max(l[2],l[1])-min(l[1],l[2]),1-max(l[2],l[1]))
    for (i in 1:(3^k-1)){
        A = rbind(A,c(min(l[2*i+1],l[2*i+2]),max(l[2*i+2],l[2*i+1])-min(l[2*i+1],l[2*i+2]),1-max(l[2*i+2],l[2*i+1])))
    }
    B<-SIMULA(A,N1)
    
    PI <- matrix(0,10,nrow=N2)
    for (i in 1:N2)
    {
        PI[i,1:10] <- ESTIMA2(B[1:(i*N1/N2)])
    }
    return(PI)
}

A = RODAR2(3,1000)
for(i in 1:20) {
    plot(type="l",y=A[i,1:5],col=i,axes=ifelse(i==1,T,F),xlab="Ordem",ylab="BIC",main="BIC para N1 = 1000 e k = 3",x=1:5)
    par(new=ifelse(i==20,F,T))
}

A = RODAR2(4,10000)
for(i in 1:20) {
    plot(type="l",y=A[i,2:6],col=i,axes=ifelse(i==1,T,F),xlab="Ordem",ylab="BIC",main="BIC para N1 = 10000 e k = 4",x=2:6)
    par(new=ifelse(i==20,F,T))
}

A = RODAR2(5,10000)
for(i in 1:20) {
    plot(type="l",y=A[i,3:7],col=i,axes=ifelse(i==1,T,F),xlab="Ordem",ylab="BIC",main="BIC para N1 = 10000 e k = 5",x=3:7)
    par(new=ifelse(i==20,F,T))
}    

A = RODAR2(6,100000)
for(i in 1:20) {
    plot(type="l",A[i,4:8],col=i,axes=ifelse(i==1,T,F),xlab="Ordem",ylab="BIC",main="BIC para N1 = 100000 e k = 6",,x=4:8)
    par(new=ifelse(i==20,F,T))
}    
