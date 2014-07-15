#É necessário instalar o pacote expm, que multiplica matrizes.
library(expm)

#função auxiliar que compara matrizes
matequal <- function(x, y)
	is.matrix(x) && is.matrix(y) && all(dim(x) == dim(y)) && all(x ==y)

#função de simulação
SIMULA <- function(A,N)
{    
    k = round(log(nrow(A),3)) #Para k muito grande o R estava dando erro de aproximação, por isso o round
    #roda os primeiros k passos, que são uniformes discretas para os valores {0,1,2}
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
    SAIDA <- SAIDA
    return (SAIDA)
}

#Função que estima a ordem, usando o BIC conforme paper http://arxiv.org/pdf/0910.0264v5.pdf
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
    #RESP = 1:10
    #return(RESP[BIC==min(BIC)])
    return(BIC)
}

RODAR <- function(k,N1,N2=20){
    l = runif(3^k*2) #A matriz usada para a simulação é gerada aleatoriamente
    A = c(min(l[1],l[2]),max(l[2],l[1])-min(l[1],l[2]),1-max(l[2],l[1]))
    for (i in 1:(3^k-1)){
        A = rbind(A,c(min(l[2*i+1],l[2*i+2]),max(l[2*i+2],l[2*i+1])-min(l[2*i+1],l[2*i+2]),1-max(l[2*i+2],l[2*i+1])))
    }
    B<-SIMULA(A,N1) #N1 é o número de passos simulados
    
    PI <- matrix(0,1,nrow=N2) #N2 é o número de cortes, para ver em qual parte do N1 é o BIC estima corretamente a ordem
    for (i in 1:N2)
    {
        PI[i] <- ESTIMA(B[1:(i*N1/N2)])
    }
    return(PI)
}

A = RODAR(3,1000)
RODAR(4,10000)
RODAR(5,10000)
RODAR(6,100000)
RODAR(7,100000)
RODAR(8,1000000)
