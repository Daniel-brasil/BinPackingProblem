#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "gurobi_c.h"

int  timeLimitACO, timeLimitMIP, NIterations, timeLimitRMP;
int  o, sr;
float QAnt, evaporation, QIncrease, pheromone, alfa, beta, nZao, QNZao, volume_Total, ZMin;
int volumeTotal, menorVolume, menorLado;
float teta;

//########################################struct probability#####################################

typedef struct {

	int id;
	float value[10];

}Probability;

typedef struct noProbability {

	Probability probability;
	struct noProbability* proximo;
	struct noProbability* permanentProximo;

}NoProbability;


Probability preencherProbability(int id) {

	Probability probability;

	probability.id = id;

	probability.value[0] = 0.0;
	probability.value[1] = 0.0;
	probability.value[2] = 0.0;
	probability.value[3] = 0.0;
	probability.value[4] = 0.0;
	probability.value[5] = 0.0;
	probability.value[6] = 0.0;
	probability.value[7] = 0.0;
	probability.value[8] = 0.0;
	probability.value[9] = 0.0;

	return probability;

}


NoProbability* insereProbability(NoProbability* topo, int id) {

	NoProbability* novo = malloc(sizeof(NoProbability));

	novo->probability = preencherProbability(id);
	novo->proximo = topo;
	novo->permanentProximo = topo;

	return novo;

}


void deleteProbability(NoProbability** topo, NoProbability* anterior, NoProbability* aDeletar) {

	if (*topo && aDeletar) {

		if (anterior == NULL || *topo == aDeletar) {

			*topo = aDeletar->proximo;

		}
		else {

			anterior->proximo = aDeletar->proximo;

		}

		aDeletar->proximo = aDeletar->permanentProximo;

	}

}

void freeProbabilities(NoProbability* topoProbability) {

	NoProbability* auxProbability;

	while (topoProbability) {

		auxProbability = topoProbability->proximo;
		free(topoProbability);
		topoProbability = auxProbability;

	}

}

//*************** Estrutura Pack ***********************
typedef struct {

	int id;
	int x;
	int y;
	int z;
	int orient;
	int kx; // projeção no eixo x
	int ky;
	int kz;

}Pack;

typedef struct noPack {

	Pack pack;
	struct noPack* proximo; //ponteiro apontando para um estrutura do tipo no

}NoPack;


//*************** Estrutura Espaço Residual ***********************
typedef struct {

	int volume;
	int x;
	int y;
	int z;
	int maxX; //tamanho maximo na eixo x
	int maxY;
	int maxZ;
	int indexBin;

}Space;


typedef struct noSpace {

	Space space;
	struct noSpace* proximoNaBin;
	struct noSpace* anteriorNaBin;
	struct noSpace* proximoGeral;
	struct noSpace* anteriorGeral;

}NoSpace;


//*************** Estrutura Bin ***********************
typedef struct {

	int idt;
	int idtColuna;
	int usedSpace;
	int qtdeItens;
	struct noPack* conteudo;
	struct noSpace* spaces;


}Bin;


typedef struct noBin {

	Bin bin;
	struct noBin* proximo; //ponteiro apontando para um estrutura do tipo noBin
	struct noBin* proximoColuna; //ponteiro apontando para um estrutura do tipo noBin

}NoBin;


//*************** Estrutura Soluction ***********************
typedef struct {

	int id;
	int value;
	float utilization;
	int time;
	struct noBin* bins;

}Soluction;


typedef struct noSoluction {

	Soluction soluction;
	struct noSoluction* proximo; //ponteiro apontando para um estrutura do tipo noSoluction

}NoSoluction;


//*************** Procedimento Pack ***********************


Pack preencherPack(int id, int x, int y, int z, int orient, int kx, int ky, int kz) {

	Pack pack;

	pack.id = id;
	pack.x = x;
	pack.y = y;
	pack.z = z;
	pack.orient = orient;
	pack.kx = kx;
	pack.ky = ky;
	pack.kz = kz;

	return pack;
}


void imprimirPack(Pack pack) {

	printf("\nid: %i,  x: %i, y: %i, z: %i, orientacao: %i, %i paralelo ao eixo x, %i paralelo ao eixo y, %i paralelo ao eixo z,\n", pack.id, pack.x, pack.y, pack.z, pack.orient, pack.kx, pack.ky, pack.kz);
}

//função operação push (empilhar)


NoPack* empilharPack(NoPack* lista, int id, int x, int y, int z, int orient, int kx, int ky, int kz) {
	NoPack* novo = malloc(sizeof(NoPack));

	if (novo) {

		novo->pack = preencherPack(id, x, y, z, orient, kx, ky, kz);//salva o dado pessoa na struct novo
		novo->proximo = lista; //salva ponteiro próximo como o último topo

		return novo;

	}
	else {

		printf("Não foi possível alocar memória");
		return NULL;
	}

}

void deletarPack(NoBin* binAtual, int idItem, int volume) {

	NoPack* auxPack, * deletar;

	auxPack = binAtual->bin.conteudo;

	if (auxPack->pack.id == idItem) { // é o primeiro item

		binAtual->bin.conteudo = auxPack->proximo;

		free(auxPack);

	}
	else {//não é o primeiro item

		while (auxPack && auxPack->proximo->pack.id != idItem) {

			auxPack = auxPack->proximo;

		}

		if (auxPack) {

			binAtual->bin.usedSpace = binAtual->bin.usedSpace - volume;

			deletar = auxPack->proximo;

			auxPack->proximo = deletar->proximo;

			free(deletar);

		}
		else {

			printf("Erro: item %d nao encontrado na Bin %d", idItem, binAtual->bin.idt);

		}


	}


}


void freeMemoryPack(NoPack* topoPack) {

	NoPack* auxPack = topoPack;
	NoPack* temp;

	while (auxPack) {

		temp = auxPack->proximo;

		free(auxPack);

		auxPack = temp;

	}


}


//*************** Procedimentos Space ***********************

Space preencherSpace(int indexBin, int x, int y, int z, int maxX, int maxY, int maxZ) {

	Space space;

	space.volume = maxX * maxY * maxZ;
	space.indexBin = indexBin;
	space.x = x;
	space.y = y;
	space.z = z;
	space.maxX = maxX; //tamanho maximo na eixo x
	space.maxY = maxY;
	space.maxZ = maxZ;

	return space;
}


//função operação push (empilhar)

int empilhaSpaceNaBin(NoBin* usedBin, NoSpace* inserir) {

	NoSpace* anterior, * aux = usedBin->bin.spaces;

	if (aux == NULL) {//lista está vazia

		usedBin->bin.spaces = inserir;
		inserir->proximoNaBin = NULL;
		inserir->anteriorNaBin = NULL;

		return 1;

	}
	else if (inserir->space.z < aux->space.z ||
		(inserir->space.z == aux->space.z && inserir->space.x < aux->space.x) ||
		(inserir->space.z == aux->space.z && inserir->space.x == aux->space.x && inserir->space.y < aux->space.y)) {
		//insere no inicio

		usedBin->bin.spaces = inserir;
		aux->anteriorNaBin = inserir;
		inserir->proximoNaBin = aux;
		inserir->anteriorNaBin = NULL;

		return 1;


	}
	else {

		//insere no meio

		anterior = aux;
		aux = aux->proximoNaBin;

		while (aux && inserir->space.z > aux->space.z) {

			anterior = aux;

			aux = aux->proximoNaBin;

		}

		while (aux && inserir->space.z == aux->space.z && inserir->space.x > aux->space.x) {

			anterior = aux;

			aux = aux->proximoNaBin;

		}

		while (aux && inserir->space.z == aux->space.z && inserir->space.x == aux->space.x && inserir->space.y > aux->space.y) {

			anterior = aux;

			aux = aux->proximoNaBin;

		}

		if (aux && inserir->space.z == aux->space.z && inserir->space.y == aux->space.y && inserir->space.x == aux->space.x) {
			//ponto duplicado não  insere

			return 0;

		}
		else if (inserir->space.z == anterior->space.z && inserir->space.y == anterior->space.y && inserir->space.x == anterior->space.x) {


			return 0;


		}
		else {

			inserir->proximoNaBin = aux;
			inserir->anteriorNaBin = anterior;

			anterior->proximoNaBin = inserir;

			if (aux) aux->anteriorNaBin = inserir;

			return 1;

		}

	}

	system("pause");

	return 1;

}

void empilhaSpaceGeral6(NoSpace** listaSpaceGeral, NoSpace* inserir) {


	NoSpace* anterior, * aux = *listaSpaceGeral;

	if (aux == NULL) {//lista está vazia

		*listaSpaceGeral = inserir;
		inserir->proximoGeral = NULL;
		inserir->anteriorGeral = NULL;

	}
	else if (inserir->space.volume < aux->space.volume) {
		//insere no inicio

		*listaSpaceGeral = inserir;
		aux->anteriorGeral = inserir;
		inserir->proximoGeral = aux;
		inserir->anteriorGeral = NULL;

	}
	else {


		anterior = aux;
		aux = aux->proximoGeral;

		while (aux && inserir->space.volume > aux->space.volume) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		inserir->proximoGeral = aux;
		inserir->anteriorGeral = anterior;

		anterior->proximoGeral = inserir;

		if (aux) aux->anteriorGeral = inserir;

	}

}


void empilhaSpaceGeral5(NoSpace** listaSpaceGeral, NoSpace* inserir) {


	NoSpace* anterior, * aux = *listaSpaceGeral;

	if (aux == NULL) {//lista está vazia

		*listaSpaceGeral = inserir;
		inserir->proximoGeral = NULL;
		inserir->anteriorGeral = NULL;

	}
	else if (inserir->space.y < aux->space.y ||
		(inserir->space.y == aux->space.y && inserir->space.x < aux->space.x) ||
		(inserir->space.y == aux->space.y && inserir->space.x == aux->space.x && inserir->space.z < aux->space.z)) {
		//insere no inicio

		*listaSpaceGeral = inserir;
		aux->anteriorGeral = inserir;
		inserir->proximoGeral = aux;
		inserir->anteriorGeral = NULL;

	}
	else {


		anterior = aux;
		aux = aux->proximoGeral;

		while (aux && inserir->space.y > aux->space.y) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.y == aux->space.y && inserir->space.x > aux->space.x) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.y == aux->space.y && inserir->space.x == aux->space.x && inserir->space.z > aux->space.z) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		inserir->proximoGeral = aux;
		inserir->anteriorGeral = anterior;

		anterior->proximoGeral = inserir;

		if (aux) aux->anteriorGeral = inserir;

	}

}

void empilhaSpaceGeral4(NoSpace** listaSpaceGeral, NoSpace* inserir) {


	NoSpace* anterior, * aux = *listaSpaceGeral;

	if (aux == NULL) {//lista está vazia

		*listaSpaceGeral = inserir;
		inserir->proximoGeral = NULL;
		inserir->anteriorGeral = NULL;

	}
	else if (inserir->space.y < aux->space.y ||
		(inserir->space.y == aux->space.y && inserir->space.z < aux->space.z) ||
		(inserir->space.y == aux->space.y && inserir->space.z == aux->space.z && inserir->space.x < aux->space.x)) {
		//insere no inicio

		*listaSpaceGeral = inserir;
		aux->anteriorGeral = inserir;
		inserir->proximoGeral = aux;
		inserir->anteriorGeral = NULL;

	}
	else {


		anterior = aux;
		aux = aux->proximoGeral;

		while (aux && inserir->space.y > aux->space.y) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.y == aux->space.y && inserir->space.z > aux->space.z) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.y == aux->space.y && inserir->space.z == aux->space.z && inserir->space.x > aux->space.x) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		inserir->proximoGeral = aux;
		inserir->anteriorGeral = anterior;

		anterior->proximoGeral = inserir;

		if (aux) aux->anteriorGeral = inserir;

	}

}

void empilhaSpaceGeral3(NoSpace** listaSpaceGeral, NoSpace* inserir) {


	NoSpace* anterior, * aux = *listaSpaceGeral;

	if (aux == NULL) {//lista está vazia

		*listaSpaceGeral = inserir;
		inserir->proximoGeral = NULL;
		inserir->anteriorGeral = NULL;

	}
	else if (inserir->space.x < aux->space.x ||
		(inserir->space.x == aux->space.x && inserir->space.z < aux->space.z) ||
		(inserir->space.x == aux->space.x && inserir->space.z == aux->space.z && inserir->space.y < aux->space.y)) {
		//insere no inicio

		*listaSpaceGeral = inserir;
		aux->anteriorGeral = inserir;
		inserir->proximoGeral = aux;
		inserir->anteriorGeral = NULL;

	}
	else {


		anterior = aux;
		aux = aux->proximoGeral;

		while (aux && inserir->space.x > aux->space.x) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.x == aux->space.x && inserir->space.z > aux->space.z) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.x == aux->space.x && inserir->space.z == aux->space.z && inserir->space.y > aux->space.y) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		inserir->proximoGeral = aux;
		inserir->anteriorGeral = anterior;

		anterior->proximoGeral = inserir;

		if (aux) aux->anteriorGeral = inserir;

	}

}

void empilhaSpaceGeral2(NoSpace** listaSpaceGeral, NoSpace* inserir) {


	NoSpace* anterior, * aux = *listaSpaceGeral;

	if (aux == NULL) {//lista está vazia

		*listaSpaceGeral = inserir;
		inserir->proximoGeral = NULL;
		inserir->anteriorGeral = NULL;

	}
	else if (inserir->space.x < aux->space.x ||
		(inserir->space.x == aux->space.x && inserir->space.y < aux->space.y) ||
		(inserir->space.x == aux->space.x && inserir->space.y == aux->space.y && inserir->space.z < aux->space.z)) {
		//insere no inicio

		*listaSpaceGeral = inserir;
		aux->anteriorGeral = inserir;
		inserir->proximoGeral = aux;
		inserir->anteriorGeral = NULL;

	}
	else {


		anterior = aux;
		aux = aux->proximoGeral;

		while (aux && inserir->space.x > aux->space.x) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.x == aux->space.x && inserir->space.y > aux->space.y) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.x == aux->space.x && inserir->space.y == aux->space.y && inserir->space.z > aux->space.z) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		inserir->proximoGeral = aux;
		inserir->anteriorGeral = anterior;

		anterior->proximoGeral = inserir;

		if (aux) aux->anteriorGeral = inserir;

	}


}

void empilhaSpaceGeral1(NoSpace** listaSpaceGeral, NoSpace* inserir) {


	NoSpace* anterior, * aux = *listaSpaceGeral;

	if (aux == NULL) {//lista está vazia

		*listaSpaceGeral = inserir;
		inserir->proximoGeral = NULL;
		inserir->anteriorGeral = NULL;

	}
	else if (inserir->space.z < aux->space.z ||
		(inserir->space.z == aux->space.z && inserir->space.y < aux->space.y) ||
		(inserir->space.z == aux->space.z && inserir->space.y == aux->space.y && inserir->space.x < aux->space.x)) {
		//insere no inicio

		*listaSpaceGeral = inserir;
		aux->anteriorGeral = inserir;
		inserir->proximoGeral = aux;
		inserir->anteriorGeral = NULL;

	}
	else {


		anterior = aux;
		aux = aux->proximoGeral;

		while (aux && inserir->space.z > aux->space.z) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.z == aux->space.z && inserir->space.y > aux->space.y) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.z == aux->space.z && inserir->space.y == aux->space.y && inserir->space.x > aux->space.x) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		inserir->proximoGeral = aux;
		inserir->anteriorGeral = anterior;

		anterior->proximoGeral = inserir;

		if (aux) aux->anteriorGeral = inserir;

	}


}


void empilhaSpaceGeral(NoSpace** listaSpaceGeral, NoSpace* inserir) {


	NoSpace* anterior, * aux = *listaSpaceGeral;

	if (aux == NULL) {//lista está vazia

		*listaSpaceGeral = inserir;
		inserir->proximoGeral = NULL;
		inserir->anteriorGeral = NULL;

	}
	else if (inserir->space.z < aux->space.z ||
		(inserir->space.z == aux->space.z && inserir->space.x < aux->space.x) ||
		(inserir->space.z == aux->space.z && inserir->space.x == aux->space.x && inserir->space.y < aux->space.y)) {
		//insere no inicio

		*listaSpaceGeral = inserir;
		aux->anteriorGeral = inserir;
		inserir->proximoGeral = aux;
		inserir->anteriorGeral = NULL;

	}
	else {


		anterior = aux;
		aux = aux->proximoGeral;

		while (aux && inserir->space.z > aux->space.z) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.z == aux->space.z && inserir->space.x > aux->space.x) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		while (aux && inserir->space.z == aux->space.z && inserir->space.x == aux->space.x && inserir->space.y > aux->space.y) {

			anterior = aux;

			aux = aux->proximoGeral;

		}

		inserir->proximoGeral = aux;
		inserir->anteriorGeral = anterior;

		anterior->proximoGeral = inserir;

		if (aux) aux->anteriorGeral = inserir;

	}


}

void empilharSpace(NoSpace** listaSpaceGeral, NoBin* usedBin, NoSpace* inserir) {

	int inserido = empilhaSpaceNaBin(usedBin, inserir);

	if (inserido == 1) {

		if (sr == 0) empilhaSpaceGeral6(listaSpaceGeral, inserir);
		else if (sr == 1) empilhaSpaceGeral(listaSpaceGeral, inserir);
		else if (sr == 2) empilhaSpaceGeral1(listaSpaceGeral, inserir);
		else if (sr == 3) empilhaSpaceGeral2(listaSpaceGeral, inserir);
		else if (sr == 4) empilhaSpaceGeral3(listaSpaceGeral, inserir);
		else if (sr == 5) empilhaSpaceGeral4(listaSpaceGeral, inserir);
		else empilhaSpaceGeral5(listaSpaceGeral, inserir);

	}
	else {

		free(inserir);

	}

}

NoSpace* criarSpace(int x, int y, int z, int indexBin, int maxX, int maxY, int maxZ) {

	NoSpace* novo = malloc(sizeof(NoSpace));

	if (novo) {

		novo->space = preencherSpace(indexBin, x, y, z, maxX, maxY, maxZ);
		novo->proximoNaBin = NULL;
		novo->proximoGeral = NULL;
		novo->anteriorNaBin = NULL;
		novo->anteriorGeral = NULL;


	}
	else {

		printf("Não foi possível alocar memória");
	}

	return novo;
}

void deletarSpace(NoSpace** topoSpaceGeral, NoBin* usedBin, NoSpace* deletar) {

	NoSpace* proximoNaBin = deletar->proximoNaBin;
	NoSpace* anteriorNaBin = deletar->anteriorNaBin;
	NoSpace* proximoGeral = deletar->proximoGeral;
	NoSpace* anteriorGeral = deletar->anteriorGeral;

	free(deletar);

	if (anteriorNaBin) {

		anteriorNaBin->proximoNaBin = proximoNaBin;

	}
	else {

		usedBin->bin.spaces = proximoNaBin;

	}

	if (proximoNaBin) {

		proximoNaBin->anteriorNaBin = anteriorNaBin;

	}


	if (anteriorGeral) {

		anteriorGeral->proximoGeral = proximoGeral;

	}
	else {

		*topoSpaceGeral = proximoGeral;

	}

	if (proximoGeral) {

		proximoGeral->anteriorGeral = anteriorGeral;

	}

}

void freeMemorySpace(NoSpace* listaSpace) {

	NoSpace* aux;

	while (listaSpace) {

		aux = listaSpace->proximoGeral;

		free(listaSpace);

		listaSpace = aux;

	}

}
//*************** Procedimentos Bin ***********************

Bin preencherBin(int id, int usedSpace, NoPack* listaPack, NoSpace* listaSpace) {

	Bin bin;

	bin.idt = id;
	bin.idtColuna = id;
	bin.usedSpace = usedSpace;
	bin.conteudo = listaPack;
	bin.spaces = listaSpace;

	return bin;
}

//função operação push (empilhar)

NoBin* empilharBin(NoBin* listaBin, NoPack* listaPack, NoSpace* listaSpace, int id, int usedSpace) {
	NoBin* novo = malloc(sizeof(NoBin));

	if (novo) {

		novo->bin = preencherBin(id, usedSpace, listaPack, listaSpace);//salva o dado pessoa na struct novo
		novo->proximo = listaBin; //salva ponteiro próximo como o último topo
		novo->proximoColuna = NULL;

		return novo;

	}
	else {

		printf("Não foi possível alocar memória");
		return NULL;
	}

}


NoBin* findBin(int idBin, NoBin* lista) { //encontra o ponteiro para um compartimento (Bin) a partir de seu index

	NoBin* usedBin = lista;

	while (usedBin && (*usedBin).bin.idt != idBin) {

		usedBin = usedBin->proximo;

	}

	return usedBin;

}

void freeMemoryBinMIP(NoBin* topoBin) {

	NoBin* auxBin = topoBin;
	NoBin* temp;

	while (auxBin) {

		freeMemoryPack(auxBin->bin.conteudo);

		temp = auxBin->proximoColuna;

		free(auxBin);

		auxBin = temp;

	}


}

void freeMemoryBin(NoBin* topoBin) {

	NoBin* auxBin = topoBin;
	NoBin* temp;

	while (auxBin) {

		freeMemoryPack(auxBin->bin.conteudo);

		temp = auxBin->proximo;

		free(auxBin);

		auxBin = temp;

	}


}

void imprimirBin(NoBin* topoBin) {

	NoBin* aux = topoBin;
	NoPack* auxPack;

	while (aux) {


		printf("\n\nBin id %d\n", aux->bin.idt);

		auxPack = aux->bin.conteudo;

		while (auxPack) {

			imprimirPack(auxPack->pack);

			auxPack = auxPack->proximo;

		}



		aux = aux->proximo;
	}

}

//******************************procedimentos Soluction*****************************

Soluction preencherSoluction(int value, float utilization, int time, int id, NoBin* listaBin) {


	Soluction soluction;

	soluction.id = id;
	soluction.value = value;
	soluction.utilization = utilization;
	soluction.time = time;
	soluction.bins = listaBin;

	return soluction;
}

NoSoluction* empilharSoluction(NoSoluction* listaSoluction, NoBin* listaBin, int id, int value, float utilization, int time) {
	NoSoluction* novo = malloc(sizeof(NoSoluction));

	if (novo) {

		novo->soluction = preencherSoluction(value, utilization, time, id, listaBin);//salva o dado pessoa na struct novo
		novo->proximo = listaSoluction; //salva ponteiro próximo como o último topo


		return novo;

	}
	else {

		printf("Não foi possível alocar memória");
		return NULL;
	}

}



void freeMemorySoluction(NoSoluction* geralSoluction) {

	freeMemoryBin(geralSoluction->soluction.bins);

	free(geralSoluction);

}

void freeMemorySoluctionAll(NoSoluction* geralSoluction) {

	NoSoluction* aux = NULL;

	while (geralSoluction) {

		geralSoluction = aux;

		geralSoluction = geralSoluction->proximo;

		freeMemorySoluction(aux);


	}


}

void freeMemorySoluctionACO(NoSoluction* geralSoluction) {

	NoSoluction* aux = NULL;

	while (geralSoluction) {

		geralSoluction = aux;

		geralSoluction = geralSoluction->proximo;

		free(aux);


	}


}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//este grupo de  funções ordena os itens de acordo com o cálculo de próximo passo e retorna um vetor com os índices ordenados

/*
int particaoIncrease(int inicio, int fim, float* nextStep, int* ordem) {
	int i = inicio;
	int j;
	int temporario;

	for (j = inicio; j < fim; j++) {

		if (nextStep[ordem[j]] <= nextStep[ordem[fim]]) {

			temporario = ordem[i];
			ordem[i] = ordem[j];
			ordem[j] = temporario;
			i = i + 1;
		}

	}

	temporario = ordem[i];
	ordem[i] = ordem[fim];
	ordem[fim] = temporario;


	return i;

}
*/

int particaoDecrease(int inicio, int fim, float* nextStep, int* ordem) {
	int i = inicio;
	int j;
	int temporario;

	for (j = inicio; j < fim; j++) {

		if (nextStep[ordem[j]] >= nextStep[ordem[fim]]) {

			temporario = ordem[i];
			ordem[i] = ordem[j];
			ordem[j] = temporario;
			i = i + 1;
		}

	}

	temporario = ordem[i];
	ordem[i] = ordem[fim];
	ordem[fim] = temporario;


	return i;

}

/*void quicksortIncrease(int inicio, int fim, float* nextStep, int* ordem) {

	if (inicio < fim) {

		int pivo = particaoIncrease(inicio, fim, nextStep, ordem);
		quicksortIncrease(inicio, pivo - 1, nextStep, ordem);
		quicksortIncrease(pivo + 1, fim, nextStep, ordem);

	}


}

*/

void quicksortDecrease(int inicio, int fim, float* nextStep, int* ordem) {

	if (inicio < fim) {

		int pivo = particaoDecrease(inicio, fim, nextStep, ordem);
		quicksortDecrease(inicio, pivo - 1, nextStep, ordem);
		quicksortDecrease(pivo + 1, fim, nextStep, ordem);

	}


}

void ordenaValores(float* nextStep, int* ordem, int quantidade) {

	int contador;

	for (contador = 0; contador < quantidade; contador++) {

		ordem[contador] = contador;

	}


	quicksortDecrease(0, quantidade - 1, nextStep, ordem);

}

//##############################################################
//importa a instância do arquivo txt e salva no tabela itens

int retornaMaior(int a, int b, int c) {

	if (a > b) {

		if (a > c) return a;
		else return c;

	}
	else {

		if (b > c) return b;
		else return c;

	}

}

int retornaMenor(int a, int b, int c) {

	if (a < b) {

		if (a < c) return a;
		else return c;

	}
	else {

		if (b < c) return b;
		else return c;

	}

}

int retornaMeio(int a, int b, int c) {

	if ((a >= b && a <= c) ||
		(a >= c && a <= b)) return a;

	else if ((b >= a && b <= c) ||
		(b >= c && b <= a)) return b;
	else return c;

}


int importaArquivoClass9(int** itens, float** heuristicas, int classe, int quantidade, int instancia, int tamanho[]) {

	volume_Total = 0.0;
	volumeTotal = 0;
	int maiorVolume = 0;
	int volume;
	int totalTipos = 0;

	char texto[30];
	char texto1[2];
	char texto2[] = "_";
	char texto3[5];
	char texto4[3];
	char texto5[] = ".txt";

	strcpy(texto, "instances/");

	sprintf(texto1, "%d", classe);
	sprintf(texto3, "%d", quantidade);
	sprintf(texto4, "%d", instancia);

	strcat(texto, texto1);
	strcat(texto, texto2);
	strcat(texto, texto3);
	strcat(texto, texto2);
	strcat(texto, texto4);
	strcat(texto, texto5);

	printf("\n\n\t++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("\n\tInstancia utilizada %s\n\n", texto);

	int contador;

	FILE* file;
	file = fopen(texto, "r");

	int tamanhoTemp[3];

	fscanf(file, "%i", &tamanhoTemp[0]);
	fscanf(file, "%i", &tamanhoTemp[1]);
	fscanf(file, "%i", &tamanhoTemp[2]);

	tamanho[0] = retornaMaior(tamanhoTemp[0], tamanhoTemp[1], tamanhoTemp[2]);
	tamanho[1] = retornaMenor(tamanhoTemp[0], tamanhoTemp[1], tamanhoTemp[2]);
	tamanho[2] = retornaMeio(tamanhoTemp[0], tamanhoTemp[1], tamanhoTemp[2]);

	menorVolume = tamanho[0] * tamanho[1] * tamanho[2];
	menorLado = retornaMenor(tamanho[0], tamanho[1], tamanho[2]);

	for (contador = 0; contador < quantidade; contador++) {

		fscanf(file, "%i", &itens[contador][0]);
		fscanf(file, "%i", &itens[contador][1]);
		fscanf(file, "%i", &itens[contador][2]);
		fscanf(file, "%i", &itens[contador][3]);

		volume = itens[contador][1] * itens[contador][2] * itens[contador][3];

		//volume
		heuristicas[0][contador] = (float)volume;
		heuristicas[1][contador] = (float)retornaMaior(itens[contador][1], itens[contador][2], itens[contador][3]);
		heuristicas[2][contador] = (float)retornaMenor(itens[contador][1], itens[contador][2], itens[contador][3]);
		heuristicas[3][contador] = (float)retornaMenor(itens[contador][1] * itens[contador][2],
			itens[contador][1] * itens[contador][3],
			itens[contador][2] * itens[contador][3]);
		heuristicas[4][contador] = (float)retornaMaior(itens[contador][1] * itens[contador][2],
			itens[contador][1] * itens[contador][3],
			itens[contador][2] * itens[contador][3]);

		volumeTotal = volumeTotal + volume;

		if (menorLado > retornaMenor(itens[contador][1], itens[contador][2], itens[contador][3]))
			menorLado = retornaMenor(itens[contador][1], itens[contador][2], itens[contador][3]);

		if (menorVolume > volume)menorVolume = volume;

		if (maiorVolume < volume) maiorVolume = volume;

		fscanf(file, "%i", &itens[contador][5]);

		if (itens[contador][5] + 1 > totalTipos) totalTipos = itens[contador][5] + 1;

	}

	fclose(file);

	if (classe > 8) {

		nZao = QNZao * ((float)maiorVolume);

	}

	return totalTipos;

}

void importaArquivo(int** itens, float** heuristicas, int classe, int quantidade, int instancia, int tamanho[]) {

	volume_Total = 0.0;
	volumeTotal = 0;
	int maiorVolume = 0;
	int volume;

	char texto[30];
	char texto1[3];
	char texto2[] = "_";
	char texto3[5];
	char texto4[3];
	char texto5[] = ".txt";

	strcpy(texto, "instances/");

	sprintf(texto1, "%d", classe);
	sprintf(texto3, "%d", quantidade);
	sprintf(texto4, "%d", instancia);

	strcat(texto, texto1);
	strcat(texto, texto2);
	strcat(texto, texto3);
	strcat(texto, texto2);
	strcat(texto, texto4);
	strcat(texto, texto5);

	printf("\n\n\t++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("\n\tInstancia utilizada %s\n\n", texto);

	int contador;

	FILE* file;
	file = fopen(texto, "r");

	menorVolume = tamanho[0] * tamanho[1] * tamanho[2];
	menorLado = retornaMenor(tamanho[0], tamanho[1], tamanho[2]);

	for (contador = 0; contador < quantidade; contador++) {

		fscanf(file, "%i", &itens[contador][0]);
		fscanf(file, "%i", &itens[contador][1]);
		fscanf(file, "%i", &itens[contador][2]);
		fscanf(file, "%i", &itens[contador][3]);

		volume = itens[contador][1] * itens[contador][2] * itens[contador][3];

		//volume
		heuristicas[0][contador] = (float)volume;
		heuristicas[1][contador] = (float)retornaMaior(itens[contador][1], itens[contador][2], itens[contador][3]);
		heuristicas[2][contador] = (float)retornaMenor(itens[contador][1], itens[contador][2], itens[contador][3]);
		heuristicas[3][contador] = (float)retornaMenor(itens[contador][1] * itens[contador][2],
			itens[contador][1] * itens[contador][3],
			itens[contador][2] * itens[contador][3]);

		volumeTotal = volumeTotal + volume;

		if (menorLado > retornaMenor(itens[contador][1], itens[contador][2], itens[contador][3]))
			menorLado = retornaMenor(itens[contador][1], itens[contador][2], itens[contador][3]);

		if (menorVolume > volume)menorVolume = volume;

		if (classe > 8) {

			if (maiorVolume < volume) maiorVolume = volume;

			volume_Total = volume_Total +
				(((float)itens[contador][1]) / 100.0 *
					((float)itens[contador][2]) / 100.0 *
					((float)itens[contador][3]) / 100.0);

		}

	}

	fclose(file);

	if (classe > 8) {

		nZao = QNZao * ((float)maiorVolume);

	}

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//indica qual a orientacao ideal de cada objeto de acordo com suas dimensoes

void firstOrientation(int** itens, int quantidade) {

	int contador;

	for (contador = 0; contador < quantidade; contador++) {

		if (itens[contador][1] < itens[contador][2]) { //x<y

			if (itens[contador][1] < itens[contador][3]) { //x<y e x<z

				if (itens[contador][2] < itens[contador][3]) { //x<y e x<z e y<z = z y x -> x z y

					itens[contador][4] = 5;

				}
				else { //x<y e x<z e y>z = y z x -> x y z

					itens[contador][4] = 1;
				}


			}
			else { //x<y e x<z e x > z : y x z -> z y x


				itens[contador][4] = 3;

			}

		}
		else { //x>y

			if (itens[contador][2] < itens[contador][3]) { //x>y e y<z

				if (itens[contador][1] < itens[contador][3]) { //x>y e y<z e x<z: z x y -> y z x

					itens[contador][4] = 4;

				}
				else { //x>y e y<z e x<z: x z y -> y x z

					itens[contador][4] = 2;

				}


			}
			else { //x>y e y>z: x y z ->

				itens[contador][4] = 6;
			}

		}

	}

}

void resolveOverlap(NoSpace* residualSpace, NoPack* item) {

	if (item->pack.x + item->pack.kx > residualSpace->space.x &&
		item->pack.y + item->pack.ky > residualSpace->space.y &&
		item->pack.z + item->pack.kz > residualSpace->space.z &&
		item->pack.x < residualSpace->space.x + residualSpace->space.maxX &&
		item->pack.y < residualSpace->space.y + residualSpace->space.maxY &&
		item->pack.z < residualSpace->space.z + residualSpace->space.maxZ) {
		//ha overlap

		if (item->pack.x - residualSpace->space.x >= item->pack.y - residualSpace->space.y &&
			item->pack.x - residualSpace->space.x >= item->pack.z - residualSpace->space.z) {

			residualSpace->space.maxX = item->pack.x - residualSpace->space.x;

		}
		else if (item->pack.y - residualSpace->space.y >= item->pack.x - residualSpace->space.x &&
			item->pack.y - residualSpace->space.y >= item->pack.z - residualSpace->space.z) {

			residualSpace->space.maxY = item->pack.y - residualSpace->space.y;

		}
		else {

			residualSpace->space.maxZ = item->pack.z - residualSpace->space.z;

		}

		residualSpace->space.volume = residualSpace->space.maxX * residualSpace->space.maxY * residualSpace->space.maxZ;

	}


}

// verifica se um item cabe em um espaço residual de acordo com uma orietanção passada

int fitSpace(NoSpace* place, int id, int orientations, int** itens) {

	switch (orientations) {

	case 1:
		if (place->space.maxX >= itens[id][1] && place->space.maxY >= itens[id][2] && place->space.maxZ >= itens[id][3]) {

			return orientations;
		}
		else {

			orientations = 0;
			return orientations;
		}
		break;
	case 2:
		if (place->space.maxX >= itens[id][2] && place->space.maxY >= itens[id][1] && place->space.maxZ >= itens[id][3]) {

			return orientations;

		}
		else {

			orientations = 0;
			return orientations;
		}
		break;
	case 3:
		if (place->space.maxX >= itens[id][3] && place->space.maxY >= itens[id][2] && place->space.maxZ >= itens[id][1]) {

			return orientations;

		}
		else {

			orientations = 0;
			return orientations;
		}

		break;
	case 4:
		if (place->space.maxX >= itens[id][2] && place->space.maxY >= itens[id][3] && place->space.maxZ >= itens[id][1]) {

			return orientations;
		}
		else {

			orientations = 0;
			return orientations;
		}
		break;
	case 5:
		if (place->space.maxX >= itens[id][1] && place->space.maxY >= itens[id][3] && place->space.maxZ >= itens[id][2]) {

			return orientations;
		}
		else {

			orientations = 0;
			return orientations;
		}
		break;
	case 6:
		if (place->space.maxX >= itens[id][3] && place->space.maxY >= itens[id][1] && place->space.maxZ >= itens[id][2]) {

			return orientations;
		}
		else {

			orientations = 0;
			return orientations;
		}
		break;
	default:
		orientations = 0;
		return orientations;
	}


}

void sequenceRotation(int sequence[], int rule) {

	int orientation;

	if (rule < 7) orientation = rule;
	else orientation = rule - 6;

	switch (orientation) {
	case 1:
		sequence[0] = 1;
		sequence[1] = 5;
		sequence[2] = 2;
		sequence[3] = 4;
		sequence[4] = 6;
		sequence[5] = 3;
		break;
	case 2:
		sequence[0] = 2;
		sequence[1] = 4;
		sequence[2] = 1;
		sequence[3] = 5;
		sequence[4] = 3;
		sequence[5] = 6;
		break;
	case 3:
		sequence[0] = 3;
		sequence[1] = 6;
		sequence[2] = 5;
		sequence[3] = 1;
		sequence[4] = 4;
		sequence[5] = 2;
		break;
	case 4:
		sequence[0] = 4;
		sequence[1] = 2;
		sequence[2] = 3;
		sequence[3] = 6;
		sequence[4] = 1;
		sequence[5] = 5;
		break;
	case 5:
		sequence[0] = 5;
		sequence[1] = 1;
		sequence[2] = 6;
		sequence[3] = 3;
		sequence[4] = 2;
		sequence[5] = 4;
		break;
	case 6:
		sequence[0] = 6;
		sequence[1] = 3;
		sequence[2] = 5;
		sequence[3] = 1;
		sequence[4] = 4;
		sequence[5] = 2;
		break;
	}

}

int packPackedMahvash(NoSpace* place, int id, int** itens) {

	int diferencas[3][3];

	int i, k, linha2 = -1, linha = -1, indice[3], menor = 10000, L = -1, W = -1, H = -1;

	diferencas[0][0] = place->space.maxX - itens[id][1];
	diferencas[0][1] = place->space.maxX - itens[id][2];
	diferencas[0][2] = place->space.maxX - itens[id][3];
	diferencas[1][0] = place->space.maxY - itens[id][1];
	diferencas[1][1] = place->space.maxY - itens[id][2];
	diferencas[1][2] = place->space.maxY - itens[id][3];
	diferencas[2][0] = place->space.maxZ - itens[id][1];
	diferencas[2][1] = place->space.maxZ - itens[id][2];
	diferencas[2][2] = place->space.maxZ - itens[id][3];

	for (i = 0; i < 3; i++) {

		for (k = 0; k < 3; k++) {

			if (diferencas[i][k] >= 0 && diferencas[i][k] < menor) {

				menor = diferencas[i][k];
				indice[i] = k;
				linha = i;


			}
		}


	}

	if (linha == -1) return 0;

	menor = 10000;

	for (i = 0; i < 3; i++) {

		for (k = 0; k < 3; k++) {

			if (diferencas[i][k] >= 0 && diferencas[i][k] < menor && i != linha && k != indice[linha]) {

				menor = diferencas[i][k];
				indice[i] = k;
				linha2 = i;


			}
		}


	}

	if (linha2 == -1) return 0;

	if (linha == 0 || linha2 == 0) L = indice[0];
	if (linha == 1 || linha2 == 1) W = indice[1];
	if (linha == 2 || linha2 == 2) H = indice[2];

	if (L == -1) {

		if ((W + H) == 1) L = 2;
		else if ((W + H) == 2) L = 1;
		else L = 0;

		if (diferencas[0][L] < 0) return 0;


	}
	else if (W == -1) {

		if ((L + H) == 1) W = 2;
		else if ((L + H) == 2) W = 1;
		else W = 0;

		if (diferencas[1][W] < 0) return 0;

	}
	else {

		if ((L + W) == 1) H = 2;
		else if ((L + W) == 2) H = 1;
		else H = 0;

		if (diferencas[2][H] < 0) return 0;

	}

	if (L == 0) {

		if (W == 1) { //x y z


			return 1;

		}
		else { //x z y

			return 5;
		}

	}
	else if (L == 1) {

		if (W == 0) { //y x z

			return 2;

		}
		else { //y z x

			return 4;

		}


	}
	else {

		if (W == 0) { //z x y

			return 6;

		}
		else { //z y x

			return 3;


		}


	}


}

int packPacked(NoSpace* place, int id, int** itens, int o) {

	int orientation = 0;
	int sequence[6];
	int j;

	sequenceRotation(sequence, o);

	for (j = 0; j < 6; j++) {

		orientation = fitSpace(place, id, sequence[j], itens);

		if (orientation) {

			return orientation;

		}

	}

	return orientation;

}

void definePointAxis(int pointAxis[], int packed, int x, int y, int z) {

	switch (packed) {
	case 1:
		pointAxis[0] = x;
		pointAxis[1] = y;
		pointAxis[2] = z;
		break;
	case 2:
		pointAxis[0] = y;
		pointAxis[1] = x;
		pointAxis[2] = z;
		break;
	case 3:
		pointAxis[0] = z;
		pointAxis[1] = y;
		pointAxis[2] = x;
		break;
	case 4:
		pointAxis[0] = y;
		pointAxis[1] = z;
		pointAxis[2] = x;
		break;
	case 5:
		pointAxis[0] = x;
		pointAxis[1] = z;
		pointAxis[2] = y;
		break;
	case 6:
		pointAxis[0] = z;
		pointAxis[1] = x;
		pointAxis[2] = y;
		break;
	}

}

void createNewSpaces(NoSpace** spaceXLeft, NoSpace** spaceXDown, NoSpace** spaceYDown, NoSpace** spaceYBack, NoSpace** spaceZLeft, NoSpace** spaceZBack, int pointAxis[], int x, int y, int z, int tamanho[], int indexBin) {


	*spaceXLeft = criarSpace(x + pointAxis[0], 0, z, indexBin, tamanho[0] - (x + pointAxis[0]), tamanho[1], tamanho[2] - z);

	*spaceXDown = criarSpace(x + pointAxis[0], y, 0, indexBin, tamanho[0] - (x + pointAxis[0]), tamanho[1] - y, tamanho[2]);

	*spaceYBack = criarSpace(0, y + pointAxis[1], z, indexBin, tamanho[0], tamanho[1] - (y + pointAxis[1]), tamanho[2] - z);

	*spaceYDown = criarSpace(x, y + pointAxis[1], 0, indexBin, tamanho[0] - x, tamanho[1] - (y + pointAxis[1]), tamanho[2]);

	*spaceZBack = criarSpace(0, y, z + pointAxis[2], indexBin, tamanho[0], tamanho[1] - y, tamanho[2] - (z + pointAxis[2]));

	*spaceZLeft = criarSpace(x, 0, z + pointAxis[2], indexBin, tamanho[0] - x, tamanho[1], tamanho[2] - (z + pointAxis[2]));

}


void updateNewSpaces(NoSpace* spaceXLeft, NoSpace* spaceXDown, NoSpace* spaceYDown, NoSpace* spaceYBack, NoSpace* spaceZLeft, NoSpace* spaceZBack, NoPack* listaPack) {

	NoPack* aux = listaPack;

	while (aux) {


		if ((*aux).pack.x <= (*spaceXLeft).space.x && (*aux).pack.y <= (*spaceYBack).space.y && (*aux).pack.z <= (*spaceZLeft).space.z) {


			if ((*aux).pack.y + (*aux).pack.ky <= (*spaceXDown).space.y) {
				//objeto está a esquerda

				//projetar no XLeft
				if ((*aux).pack.z <= (*spaceXLeft).space.z && (*aux).pack.z + (*aux).pack.kz > (*spaceXLeft).space.z
					&& (*aux).pack.x <= (*spaceXLeft).space.x && (*aux).pack.x + (*aux).pack.kx > (*spaceXLeft).space.x
					&& (*aux).pack.y + (*aux).pack.ky > (*spaceXLeft).space.y) {

					(*spaceXLeft).space.maxY = (*spaceXLeft).space.maxY - ((*aux).pack.y + (*aux).pack.ky - (*spaceXLeft).space.y);
					(*spaceXLeft).space.y = (*aux).pack.y + (*aux).pack.ky;

				}

				//projetar no ZLeft
				if ((*aux).pack.z <= (*spaceZLeft).space.z && (*aux).pack.z + (*aux).pack.kz > (*spaceZLeft).space.z
					&& (*aux).pack.x <= (*spaceZLeft).space.x && (*aux).pack.x + (*aux).pack.kx > (*spaceZLeft).space.x
					&& (*aux).pack.y + (*aux).pack.ky > (*spaceZLeft).space.y) {

					(*spaceZLeft).space.maxY = (*spaceZLeft).space.maxY - ((*aux).pack.y + (*aux).pack.ky - (*spaceZLeft).space.y);
					(*spaceZLeft).space.y = (*aux).pack.y + (*aux).pack.ky;

				}


			}
			else if ((*aux).pack.x + (*aux).pack.kx <= (*spaceYDown).space.x) {
				//objeto está atrás

				//projetar no YBack
				if ((*aux).pack.z <= (*spaceYBack).space.z && (*aux).pack.z + (*aux).pack.kz > (*spaceYBack).space.z
					&& (*aux).pack.y <= (*spaceYBack).space.y && (*aux).pack.y + (*aux).pack.ky > (*spaceYBack).space.y
					&& (*aux).pack.x + (*aux).pack.kx > (*spaceYBack).space.x) {

					(*spaceYBack).space.maxX = (*spaceYBack).space.maxX - ((*aux).pack.x + (*aux).pack.kx - (*spaceYBack).space.x);
					(*spaceYBack).space.x = (*aux).pack.x + (*aux).pack.kx;

				}


				//projetar no ZBack
				if ((*aux).pack.z <= (*spaceZBack).space.z && (*aux).pack.z + (*aux).pack.kz > (*spaceZBack).space.z
					&& (*aux).pack.y <= (*spaceZBack).space.y && (*aux).pack.y + (*aux).pack.ky > (*spaceZBack).space.y
					&& (*aux).pack.x + (*aux).pack.kx > (*spaceZBack).space.x) {

					(*spaceZBack).space.maxX = (*spaceZBack).space.maxX - ((*aux).pack.x + (*aux).pack.kx - (*spaceZBack).space.x);
					(*spaceZBack).space.x = (*aux).pack.x + (*aux).pack.kx;

				}


			}
			else if ((*aux).pack.z + (*aux).pack.kz <= (*spaceXLeft).space.z) {
				//objeto está abaixo

				//projetar no XDown
				if ((*aux).pack.x <= (*spaceXDown).space.x && (*aux).pack.x + (*aux).pack.kx > (*spaceXDown).space.x
					&& (*aux).pack.y <= (*spaceXDown).space.y && (*aux).pack.y + (*aux).pack.ky > (*spaceXDown).space.y
					&& (*aux).pack.z + (*aux).pack.kz > (*spaceXDown).space.z) {

					(*spaceXDown).space.maxZ = (*spaceXDown).space.maxZ - ((*aux).pack.z + (*aux).pack.kz - (*spaceXDown).space.z);
					(*spaceXDown).space.z = (*aux).pack.z + (*aux).pack.kz;

				}

				//projetar o YDown
				if ((*aux).pack.x <= (*spaceYDown).space.x && (*aux).pack.x + (*aux).pack.kx > (*spaceYDown).space.x
					&& (*aux).pack.y <= (*spaceYDown).space.y && (*aux).pack.y + (*aux).pack.ky > (*spaceYDown).space.y
					&& (*aux).pack.z + (*aux).pack.kz > (*spaceYDown).space.z) {

					(*spaceYDown).space.maxZ = (*spaceYDown).space.maxZ - ((*aux).pack.z + (*aux).pack.kz - (*spaceYDown).space.z);
					(*spaceYDown).space.z = (*aux).pack.z + (*aux).pack.kz;

				}


			}

		}

		//system("pause");

		aux = aux->proximo;
	}

	aux = listaPack;

	while (aux) {

		resolveOverlap(spaceXLeft, aux);

		resolveOverlap(spaceXDown, aux);

		resolveOverlap(spaceYDown, aux);

		resolveOverlap(spaceYBack, aux);

		resolveOverlap(spaceZLeft, aux);

		resolveOverlap(spaceZBack, aux);

		aux = aux->proximo;

	}

	spaceXLeft->space.volume = spaceXLeft->space.maxX * spaceXLeft->space.maxY * spaceXLeft->space.maxZ;
	spaceXDown->space.volume = spaceXDown->space.maxX * spaceXDown->space.maxY * spaceXDown->space.maxZ;
	spaceYDown->space.volume = spaceYDown->space.maxX * spaceYDown->space.maxY * spaceYDown->space.maxZ;
	spaceYBack->space.volume = spaceYBack->space.maxX * spaceYBack->space.maxY * spaceYBack->space.maxZ;
	spaceZLeft->space.volume = spaceZLeft->space.maxX * spaceZLeft->space.maxY * spaceZLeft->space.maxZ;
	spaceZBack->space.volume = spaceZBack->space.maxX * spaceZBack->space.maxY * spaceZBack->space.maxZ;

}

void reclassificaSpace(NoSpace** spaceGeral, NoSpace* space) {


	NoSpace* atual = space, * anterior = space->anteriorGeral;

	if (anterior && atual->space.volume < anterior->space.volume) {

		NoSpace* anteriorAoAnterior, * aposOProximo;

		while (anterior && atual->space.volume < anterior->space.volume) {

			anteriorAoAnterior = anterior->anteriorGeral;

			aposOProximo = atual->proximoGeral;

			if (aposOProximo) aposOProximo->anteriorGeral = anterior;
			anterior->proximoGeral = aposOProximo;

			anterior->anteriorGeral = atual;
			atual->proximoGeral = anterior;

			if (anteriorAoAnterior) {

				anteriorAoAnterior->proximoGeral = atual;

			}
			else {

				*spaceGeral = atual;

			}

			atual->anteriorGeral = anteriorAoAnterior;

			anterior = atual->anteriorGeral;

		}

	}

}

void updateSpaces(NoSpace** spaceGeral, NoBin* usedBin, int x, int y, int z, int pointAxis[3]) {

	NoPack* item = empilharPack(NULL, 0, x, y, z, 1, pointAxis[0], pointAxis[1], pointAxis[2]);

	NoSpace* deletar, * aux = usedBin->bin.spaces;

	while (aux) {

		resolveOverlap(aux, item);

		if (aux->space.volume < menorVolume || retornaMenor(aux->space.maxX, aux->space.maxY, aux->space.maxZ) < menorLado) {

			deletar = aux;

			aux = aux->proximoNaBin;

			deletarSpace(spaceGeral, usedBin, deletar);

		}
		else {


			if (sr == 0)  reclassificaSpace(spaceGeral, aux);

			aux = aux->proximoNaBin;



		}


	}


	free(item);

}


void deletRepeatedSpaces(NoSpace** spaceXLeft, NoSpace** spaceXDown, NoSpace** spaceYDown, NoSpace** spaceYBack, NoSpace** spaceZLeft, NoSpace** spaceZBack) {

	NoSpace* aux = *spaceXLeft;

	NoSpace* aux1 = *spaceXDown;

	if (aux->space.y == aux1->space.y) {

		free(aux);

		*spaceXLeft = NULL;

	}
	else if (aux->space.z == aux1->space.z) {

		free(aux1);

		*spaceXDown = NULL;

	}

	aux = *spaceYDown;

	aux1 = *spaceYBack;

	if (aux->space.x == aux1->space.x) {

		free(aux1);

		*spaceYBack = NULL;

	}
	else if (aux->space.z == aux1->space.z) {

		free(aux);

		*spaceYDown = NULL;

	}


	aux = *spaceZLeft;

	aux1 = *spaceZBack;

	if (aux->space.y == aux1->space.y) {

		free(aux);

		*spaceZLeft = NULL;

	}
	else if (aux->space.x == aux1->space.x) {

		free(aux1);

		*spaceZBack = NULL;

	}

	//exclui espaço com volume 0

	aux = *spaceXLeft;

	if (aux && (aux->space.volume <= menorVolume || retornaMenor(aux->space.maxX, aux->space.maxY, aux->space.maxZ) < menorLado)) {

		free(aux);
		*spaceXLeft = NULL;

	}

	aux = *spaceXDown;

	if (aux && (aux->space.volume <= menorVolume || retornaMenor(aux->space.maxX, aux->space.maxY, aux->space.maxZ) < menorLado)) {

		free(aux);
		*spaceXDown = NULL;

	}

	aux = *spaceYDown;

	if (aux && (aux->space.volume <= menorVolume || retornaMenor(aux->space.maxX, aux->space.maxY, aux->space.maxZ) < menorLado)) {

		free(aux);
		*spaceYDown = NULL;

	}

	aux = *spaceYBack;

	if (aux && (aux->space.volume <= menorVolume || retornaMenor(aux->space.maxX, aux->space.maxY, aux->space.maxZ) < menorLado)) {

		free(aux);
		*spaceYBack = NULL;

	}

	aux = *spaceZLeft;

	if (aux && (aux->space.volume <= menorVolume || retornaMenor(aux->space.maxX, aux->space.maxY, aux->space.maxZ) < menorLado)) {

		free(aux);
		*spaceZLeft = NULL;

	}

	aux = *spaceZBack;

	if (aux && (aux->space.volume <= menorVolume || retornaMenor(aux->space.maxX, aux->space.maxY, aux->space.maxZ) < menorLado)) {

		free(aux);
		*spaceZBack = NULL;

	}

}


NoBin* binPack(int** itens, int* ordem, int quantidade, int tamanho[], double* PiValue, int classe) {

	int i, x, y, z, id, packed = 0;
	int idBin;
	int pointAxis[3];
	NoBin* topoBin, * usedBin;
	NoSpace* aux, * topoSpaceGeral;
	NoPack* topoPack;
	NoSpace* spaceXLeft, * spaceXDown, * spaceYDown, * spaceYBack, * spaceZLeft, * spaceZBack, * novoSpace;
	double valorDual;

	topoBin = NULL;
	topoSpaceGeral = NULL;
	topoPack = NULL;
	idBin = 0;

	topoBin = empilharBin(topoBin, topoPack, NULL, idBin, 0);
	topoBin->bin.qtdeItens = 0;

	novoSpace = criarSpace(0, 0, 0, idBin, tamanho[0], tamanho[1], tamanho[2]);
	if (novoSpace) empilharSpace(&topoSpaceGeral, topoBin, novoSpace);

	idBin++;

	for (i = 0; i < quantidade; i++) {

		id = ordem[i];

		packed = 0; //informa orientação do empacotamento, se não empacotado, informa 0
		aux = topoSpaceGeral;

		while (aux && packed == 0) {

			if (aux->space.volume >= itens[id][1] * itens[id][2] * itens[id][3]) {

				if (o == -1) packed = packPackedMahvash(aux, id, itens);
				else packed = packPacked(aux, id, itens, o + itens[id][4]);

				if (packed) {

					spaceXLeft = NULL;
					spaceXDown = NULL;
					spaceYDown = NULL;
					spaceYBack = NULL;
					spaceZLeft = NULL;
					spaceZBack = NULL;

					definePointAxis(pointAxis, packed, itens[id][1], itens[id][2], itens[id][3]);

					usedBin = findBin((*aux).space.indexBin, topoBin);

					x = (*aux).space.x;
					y = (*aux).space.y;
					z = (*aux).space.z;

					createNewSpaces(&spaceXLeft, &spaceXDown, &spaceYDown, &spaceYBack, &spaceZLeft, &spaceZBack, pointAxis, x, y, z, tamanho, aux->space.indexBin);

					updateNewSpaces(spaceXLeft, spaceXDown, spaceYDown, spaceYBack, spaceZLeft, spaceZBack, (*usedBin).bin.conteudo);

					deletRepeatedSpaces(&spaceXLeft, &spaceXDown, &spaceYDown, &spaceYBack, &spaceZLeft, &spaceZBack);

					deletarSpace(&topoSpaceGeral, usedBin, aux);

					updateSpaces(&topoSpaceGeral, usedBin, x, y, z, pointAxis);

					if (spaceZLeft) empilharSpace(&topoSpaceGeral, usedBin, spaceZLeft);
					if (spaceZBack) empilharSpace(&topoSpaceGeral, usedBin, spaceZBack);
					if (spaceXLeft) empilharSpace(&topoSpaceGeral, usedBin, spaceXLeft);
					if (spaceYBack) empilharSpace(&topoSpaceGeral, usedBin, spaceYBack);
					if (spaceXDown) empilharSpace(&topoSpaceGeral, usedBin, spaceXDown);
					if (spaceYDown) empilharSpace(&topoSpaceGeral, usedBin, spaceYDown);

					topoPack = empilharPack((*usedBin).bin.conteudo, id, x, y, z, packed, pointAxis[0], pointAxis[1], pointAxis[2]);

					(*usedBin).bin.conteudo = topoPack;
					(*usedBin).bin.usedSpace = (*usedBin).bin.usedSpace + (itens[id][1] * itens[id][2] * itens[id][3]);
					usedBin->bin.qtdeItens = usedBin->bin.qtdeItens + 1;

				}

			}

			aux = aux->proximoGeral;

		}

		if (packed == 0) {

			if (classe == 9) {

				valorDual = PiValue[itens[id][5]];

			}
			else {

				valorDual = PiValue[id];

			}

			if (valorDual > 0.0) {

				topoPack = NULL;

				topoBin = empilharBin(topoBin, topoPack, NULL, idBin, 0);
				topoBin->bin.qtdeItens = 0;

				novoSpace = criarSpace(0, 0, 0, idBin, tamanho[0], tamanho[1], tamanho[2]);

				if (o == -1) packed = packPackedMahvash(novoSpace, id, itens);
				else packed = packPacked(novoSpace, id, itens, o + itens[id][4]);

				free(novoSpace);
				novoSpace = NULL;

				if (packed) {

					definePointAxis(pointAxis, packed, itens[id][1], itens[id][2], itens[id][3]);

					novoSpace = criarSpace(0, 0, pointAxis[2], idBin, tamanho[0], tamanho[1], tamanho[2] - pointAxis[2]);
					if (novoSpace) empilharSpace(&topoSpaceGeral, topoBin, novoSpace);

					novoSpace = criarSpace(pointAxis[0], 0, 0, idBin, tamanho[0] - pointAxis[0], tamanho[1], tamanho[2]);
					if (novoSpace) empilharSpace(&topoSpaceGeral, topoBin, novoSpace);

					novoSpace = criarSpace(0, pointAxis[1], 0, idBin, tamanho[0], tamanho[1] - pointAxis[1], tamanho[2]);
					if (novoSpace) empilharSpace(&topoSpaceGeral, topoBin, novoSpace);

					topoPack = empilharPack((*topoBin).bin.conteudo, id, 0, 0, 0, packed, pointAxis[0], pointAxis[1], pointAxis[2]);

					(*topoBin).bin.conteudo = topoPack;
					(*topoBin).bin.usedSpace = (*topoBin).bin.usedSpace + (itens[id][1] * itens[id][2] * itens[id][3]);
					topoBin->bin.qtdeItens = (*topoBin).bin.qtdeItens + 1;

					idBin++;

				}
				else {

					printf("\nItem %d nao carregado\n", id);
					system("pause");

				}

			}


		}

	}

	freeMemorySpace(topoSpaceGeral);

	return topoBin;

}


//define o tamanho da Bin de acordo com a classe do problema
void tamanhoInstancia(int classe, int tamanho[]) {


	if (classe == 1 || classe == 2 || classe == 3 || classe == 4 || classe == 5 || classe == 8) {

		tamanho[0] = 100;
		tamanho[1] = 100;
		tamanho[2] = 100;
		nZao = QNZao * 1000000.0;

	}
	else if (classe == 6) {

		tamanho[0] = 10;
		tamanho[1] = 10;
		tamanho[2] = 10;
		nZao = QNZao * 1000.0;

	}
	else if (classe == 7) {

		tamanho[0] = 40;
		tamanho[1] = 40;
		tamanho[2] = 40;
		nZao = QNZao * 42875.0;

	}
	else if (classe == 10) {

		tamanho[0] = 609;
		tamanho[1] = 243;
		tamanho[2] = 243;

	}

}


int validSoluctionClass9(NoBin* listaBin, int tam[], int solucao, int totalTipos, int* qtdTipos, int** itens) {

	int erro = 0, totalBin = 0, tipo, k;
	int* used = (int*)calloc(totalTipos, sizeof(int));
	NoBin* auxBin = listaBin;
	NoPack* auxPack;
	NoPack* auxPack1;

	//int volumeTotal;

	while (auxBin) {// percorre todas as bins utilizadas

		totalBin = totalBin + auxBin->bin.idt;

		//volumeTotal = 0;

		//printf("\nBin %d\n", auxBin->bin.idt);

		auxPack = (*auxBin).bin.conteudo;

		while (auxPack) { // percorre todos os elementos da bin


			//volumeTotal = volumeTotal + (auxPack->pack.kx * auxPack->pack.ky * auxPack->pack.kz);

			tipo = itens[auxPack->pack.id][5];

			used[tipo] = used[tipo] + auxBin->bin.idt;

			if ((*auxPack).pack.x + (*auxPack).pack.kx > tam[0] || (*auxPack).pack.y + (*auxPack).pack.ky > tam[1] || (*auxPack).pack.z + (*auxPack).pack.kz > tam[2]) {

				printf("\n\tSolucao invalida! Objeto ultrapassa os limites da bin!\n");
				imprimirPack(auxPack->pack);

				erro++;

			}


			auxPack1 = auxPack->proximo;

			while (auxPack1) { // percorre todos os demais elementos da bin

				if ((*auxPack).pack.x + (*auxPack).pack.kx > (*auxPack1).pack.x &&
					(*auxPack1).pack.x + (*auxPack1).pack.kx > (*auxPack).pack.x &&
					(*auxPack).pack.y + (*auxPack).pack.ky > (*auxPack1).pack.y &&
					(*auxPack1).pack.y + (*auxPack1).pack.ky > (*auxPack).pack.y &&
					(*auxPack).pack.z + (*auxPack).pack.kz > (*auxPack1).pack.z &&
					(*auxPack1).pack.z + (*auxPack1).pack.kz > (*auxPack).pack.z) {

					printf("\n\tSolucao invalida! Objetos com sobreposicao!\n");
					imprimirPack(auxPack->pack);
					imprimirPack(auxPack1->pack);

					erro++;

				}


				auxPack1 = auxPack1->proximo;

			}


			auxPack = auxPack->proximo;
		}

		//printf("\nVolume Total carregado %d\n", volumeTotal);

		auxBin = auxBin->proximo;
	}

	for (k = 0; k < totalTipos; k++) {

		if (qtdTipos[k] > used[k]) {

			printf("\n\tItens tipo %d empacotados %d! Itens totais %d!\n", k, used[k], qtdTipos[k]);

			erro++;


		}

	}

	if (totalBin != solucao) {

		printf("\n\tBins utilizadas %d! Solucao %d!\n", totalBin, solucao);

		erro++;

	}

	free(used);

	return erro;


}


//valida solução final para verificar todos os itens estão no limite da Bin e se não há overleap entre nenhum item
int validSoluction(NoBin* listaBin, int tam[], int quantidade, int solucao) {

	int erro = 0, totalItens = 0, totalBin = 0;
	int* used = (int*)calloc(quantidade, sizeof(int));
	NoBin* auxBin = listaBin;
	NoPack* auxPack;
	NoPack* auxPack1;

	//int volumeTotal;

	while (auxBin) {// percorre todas as bins utilizadas

		totalBin++;

		//volumeTotal = 0;

		//printf("\nBin %d\n", auxBin->bin.idt);

		auxPack = (*auxBin).bin.conteudo;

		while (auxPack) { // percorre todos os elementos da bin


			//volumeTotal = volumeTotal + (auxPack->pack.kx * auxPack->pack.ky * auxPack->pack.kz);

			if (used[auxPack->pack.id] == 0) {

				used[auxPack->pack.id] = 1;

			}
			else {

				printf("\n\tItem %d duplicado!\n", auxPack->pack.id);

				erro++;


			}


			totalItens++;

			if ((*auxPack).pack.x + (*auxPack).pack.kx > tam[0] || (*auxPack).pack.y + (*auxPack).pack.ky > tam[1] || (*auxPack).pack.z + (*auxPack).pack.kz > tam[2]) {

				printf("\n\tSolucao invalida! Objeto ultrapassa os limites da bin!\n");
				imprimirPack(auxPack->pack);

				erro++;

			}


			auxPack1 = auxPack->proximo;

			while (auxPack1) { // percorre todos os demais elementos da bin

				if ((*auxPack).pack.x + (*auxPack).pack.kx > (*auxPack1).pack.x &&
					(*auxPack1).pack.x + (*auxPack1).pack.kx > (*auxPack).pack.x &&
					(*auxPack).pack.y + (*auxPack).pack.ky > (*auxPack1).pack.y &&
					(*auxPack1).pack.y + (*auxPack1).pack.ky > (*auxPack).pack.y &&
					(*auxPack).pack.z + (*auxPack).pack.kz > (*auxPack1).pack.z &&
					(*auxPack1).pack.z + (*auxPack1).pack.kz > (*auxPack).pack.z) {

					printf("\n\tSolucao invalida! Objetos com sobreposicao!\n");
					imprimirPack(auxPack->pack);
					imprimirPack(auxPack1->pack);

					erro++;

				}


				auxPack1 = auxPack1->proximo;

			}


			auxPack = auxPack->proximo;
		}

		//printf("\nVolume Total carregado %d\n", volumeTotal);

		auxBin = auxBin->proximo;
	}

	if (totalItens != quantidade) {

		printf("\n\tItens empacotados %d! Itens totais %d!\n", totalItens, quantidade);

		erro++;

	}

	if (totalBin != solucao) {

		printf("\n\tBins utilizadas %d! Solucao %d!\n", totalBin, solucao);

		erro++;

	}

	free(used);

	return erro;

}


//salva o arquivo com os resultados
void saveFile(NoBin* listaBin, int classe, int quantidade, int instancia, int objetivo, int totalTime, unsigned int semente) {

	FILE* pont_arqu;
	NoBin* auxBin = listaBin;
	NoPack* auxPack;

	float wasted = 0.0;

	if (classe == 10) {

		wasted = ((float)objetivo) * 36.24556364 - volume_Total;

	}

	pont_arqu = fopen("BPP-0.1.6.csv", "a");

	fprintf(pont_arqu, "%d;%d;%d;%d;%d;%d;%2f;\n", classe, quantidade, instancia, semente, objetivo, totalTime, wasted);

	/*while (auxBin) {// percorre todas as bins utilizadas

		auxPack = (*auxBin).bin.conteudo;

		fprintf(pont_arqu, "[");

		while (auxPack) { // percorre todos os elementos da bin

			fprintf(pont_arqu, "[%d, %d, %d, %d, %d]", (*auxPack).pack.id, (*auxPack).pack.x, (*auxPack).pack.y, (*auxPack).pack.z, (*auxPack).pack.orient);

			if (auxPack->proximo) {

				fprintf(pont_arqu, ",");

			}

			auxPack = auxPack->proximo;
		}

		if (auxBin->proximo) {

			fprintf(pont_arqu, "],");

		}
		else {

			fprintf(pont_arqu, "];\n");
		}


		auxBin = auxBin->proximo;
	}
*/
	fclose(pont_arqu);

}


/*void defineClassificationType(int **itens, int quantidade, int *nextStep){

	int i;

	for(i=0; i<quantidade; i++){

		if(classification==0){//volume

			nextStep[i] = itens[i][1]*itens[i][2]*itens[i][3];

		}else if(classification==1){//maior lado

			nextStep[i] = retornaMaior(itens[i][1], itens[i][2], itens[i][3]);

		}else if(classification==2){//menor lado

			nextStep[i] = retornaMenor(itens[i][1], itens[i][2], itens[i][3]);

		}else if(classification==3){//maior area

			nextStep[i] = retornaMaior(itens[i][1]*itens[i][2], itens[i][1]*itens[i][3], itens[i][2]*itens[i][3]);

		}else if(classification==4){//menor area

			nextStep[i] = retornaMenor(itens[i][1]*itens[i][2], itens[i][1]*itens[i][3], itens[i][2]*itens[i][3]);

		}else if(classification==5){//menor diferenca de lado

			nextStep[i] = retornaMenor( abs(itens[i][1]-itens[i][2]),
										abs(itens[i][1]-itens[i][3]),
										abs(itens[i][2]-itens[i][3]));

		}else if(classification==6){//maior diferenca de lado

			nextStep[i] = retornaMaior( abs(itens[i][1]-itens[i][2]),
										abs(itens[i][1]-itens[i][3]),
										abs(itens[i][2]-itens[i][3]));

		}else if(classification==7){//maior diferenca de area

			nextStep[i] = retornaMaior( abs(itens[i][1]*itens[i][2]-itens[i][1]*itens[i][3]),
										abs(itens[i][1]*itens[i][2]-itens[i][2]*itens[i][3]),
										abs(itens[i][1]*itens[i][3]-itens[i][2]*itens[i][3]));

		}else{//menor diferenca de area

			nextStep[i] = retornaMenor( abs(itens[i][1]*itens[i][2]-itens[i][1]*itens[i][3]),
										abs(itens[i][1]*itens[i][2]-itens[i][2]*itens[i][3]),
										abs(itens[i][1]*itens[i][3]-itens[i][2]*itens[i][3]));
		}

	}

}
*/

void printOrdem(int** itens, int* nextStep, int* ordem, int quantidade) {

	int i;

	for (i = 0; i < quantidade; i++) {

		printf("\n%d\t%d\t%d\t%d\t%d", ordem[i], itens[ordem[i]][1], itens[ordem[i]][2], itens[ordem[i]][3], nextStep[ordem[i]]);

	}

}


/*void defineOrientationRuleAndSpaceRule(int classe){

	if(classe<9){

		orientationRule = 1;
		spaceRule = 2;

	}else{

		orientationRule = 6;
		spaceRule = 7;

	}

}
*/

void preencheCombinacoes(int** combinacoes, int quantidadeCombinacoes) {

	if (quantidadeCombinacoes > 0) { combinacoes[0][0] = 0;	combinacoes[0][1] = -1;	combinacoes[0][2] = 0; }

	if (quantidadeCombinacoes > 1) { combinacoes[1][0] = 0;	combinacoes[1][1] = 0;	combinacoes[1][2] = 0; }

	if (quantidadeCombinacoes > 2) { combinacoes[2][0] = 1;	combinacoes[2][1] = 5;	combinacoes[2][2] = 3; }

	if (quantidadeCombinacoes > 3) { combinacoes[3][0] = 1;	combinacoes[3][1] = 5;	combinacoes[3][2] = 4; }

	if (quantidadeCombinacoes > 4) { combinacoes[4][0] = 2;	combinacoes[4][1] = -1;	combinacoes[4][2] = 0; }

	if (quantidadeCombinacoes > 5) { combinacoes[5][0] = 3;	combinacoes[5][1] = -1;	combinacoes[5][2] = 0; }

	if (quantidadeCombinacoes > 6) { combinacoes[6][0] = 3;	combinacoes[6][1] = 0;	combinacoes[6][2] = 0; }

	if (quantidadeCombinacoes > 7) { combinacoes[7][0] = 0;	combinacoes[7][1] = -1;	combinacoes[7][2] = 3; }

	if (quantidadeCombinacoes > 8) { combinacoes[8][0] = 0;	combinacoes[8][1] = 2;	combinacoes[8][2] = 5; }

	if (quantidadeCombinacoes > 9) { combinacoes[9][0] = 4;	combinacoes[9][1] = 4;	combinacoes[9][2] = 4; }

}



float calcSmallerUtilization(NoBin* topoBin, int tamanho[]) {

	NoBin* aux = topoBin;
	int auxUtilization = aux->bin.usedSpace;
	float utilization;

	while (aux) {

		if (aux->bin.usedSpace < auxUtilization) {

			auxUtilization = aux->bin.usedSpace;
		}

		aux = aux->proximo;

	}

	utilization = (float)auxUtilization / ((float)(tamanho[0] * tamanho[1] * tamanho[2]));

	return utilization;

}

//############################ Ant Colony Optmization ####################################


int sorteiaProbabilidadeCombinacao(float* probabilityCombination, int quantidadeCombinacoes) {

	float luck;
	int j = 0;
	float acumulated = probabilityCombination[0];

	luck = ((float)rand() / (float)RAND_MAX);

	while (j < quantidadeCombinacoes - 1 && luck < acumulated + probabilityCombination[j + 1]) {

		j++;
		acumulated = acumulated + probabilityCombination[j];

	}

	return j;

}


void calculaProbabilityCombination(float* probabilityCombination, float totalProbabilityCombination, int quantidadeCombinacoes) {


	int k;

	for (k = 0; k < quantidadeCombinacoes; k++) {

		probabilityCombination[k] = probabilityCombination[k] / totalProbabilityCombination;

	}


}

NoProbability* montaPilhaProbabilidade(int quantidade, float** pheromones, float** heuristicas, int quantidadeCombinacoes, int quantidadeOrdens) {

	NoProbability* aux = NULL;
	int i, k;

	for (i = quantidade - 1; i >= 0; i--) {

		aux = insereProbability(aux, i);

		for (k = 0; k < quantidadeCombinacoes; k++) {

			pheromones[k][i] = pheromone;

		}

		for (k = 0; k < quantidadeOrdens; k++) {

			if (k == 0) { heuristicas[k][i] = (float)pow((double)(heuristicas[k][i] / nZao), (double)beta); }
			else { heuristicas[k][i] = (float)pow((double)(heuristicas[k][i]), (double)beta); }

		}

	}

	return aux;

}

//calcula o denominador. Calculado a cada itera��o
void calcDenominator(int** combinacoes, float** heuristicas, float** pheromones, int quantidade, NoProbability* topo, float* denominator, int quantidadeCombinacoes) {

	NoProbability* aux = topo;
	int j, k;


	for (k = 0; k < quantidadeCombinacoes; k++) {

		denominator[k] = 0.0;

	}

	for (j = 0; j < quantidade; j++) {

		for (k = 0; k < quantidadeCombinacoes; k++) {

			aux->probability.value[k] = (float)pow((double)pheromones[k][j], (double)alfa) * heuristicas[combinacoes[k][0]][j];

			denominator[k] = denominator[k] + aux->probability.value[k];

		}

		aux = aux->permanentProximo;

	}


}

//calcula a probabilidade, deve ser calculado para cada formiga
void calcProbability(int quantidade, float denominator[7], NoProbability* topo, int* ordem, int k) {

	int j;
	float luck, acumulated;
	NoProbability* aux, * anterior, * novotopo = topo;
	float newdenominator = denominator[k];

	for (j = 0; j < quantidade; j++) {

		luck = ((float)rand() / (float)RAND_MAX);

		aux = novotopo;

		anterior = NULL;

		acumulated = aux->probability.value[k] / newdenominator;

		while (aux->proximo && acumulated < luck) {

			anterior = aux;

			aux = aux->proximo;

			acumulated = acumulated + (aux->probability.value[k] / newdenominator);

		}

		ordem[j] = aux->probability.id;

		newdenominator = newdenominator - aux->probability.value[k];

		deleteProbability(&novotopo, anterior, aux);

	}

}

void copyVector(int* ordem, int* bestOrder, int quantidade) {

	int i;

	for (i = 0; i < quantidade; i++) {

		bestOrder[i] = ordem[i];

	}
}

void updatePheromone(int* bestOrder, float* pheromones, int quantidade, float increase, float decrease) {

	int i;

	for (i = 0; i < quantidade; i++) {


		pheromones[bestOrder[i]] = (pheromones[bestOrder[i]] * (((float)1.0) - evaporation)) + pheromone - (((float)i) * decrease) + ((float)increase) * QIncrease;

	}

}

int updateProbabilityCombination(float* probabilityCombination, float totalProbabilityCombination, int bestCombination, float increaseValue, int quantidadeCombinacoes) {

	int i;

	float newTotalProbabilityCombination = 0.0;

	for (i = 0; i < quantidadeCombinacoes; i++) {

		probabilityCombination[i] = probabilityCombination[i] * totalProbabilityCombination * ((float)1.0 - evaporation);

		newTotalProbabilityCombination = newTotalProbabilityCombination + probabilityCombination[i];

	}

	probabilityCombination[bestCombination] = probabilityCombination[bestCombination] + (increaseValue / totalProbabilityCombination);

	newTotalProbabilityCombination = newTotalProbabilityCombination + (increaseValue / totalProbabilityCombination);

	for (i = 0; i < quantidadeCombinacoes; i++) {

		probabilityCombination[i] = probabilityCombination[i] / newTotalProbabilityCombination;

	}

	return newTotalProbabilityCombination;

}

//############################# par�metros ###################################

void printParams() {

	printf("\n\tValores dos parametros");
	printf("\n\tQAnt\t\t\t %f", QAnt);
	printf("\n\tAlpha value\t\t %f", alfa);
	printf("\n\tBeta value\t\t %f", beta);
	printf("\n\tEvaporation proportion\t %f", evaporation);
	printf("\n\tAddition proportion\t %f", QIncrease);
	printf("\n\tPheromone value\t\t %f", pheromone);
	printf("\n\tQNZao\t\t\t %f", QNZao);
	printf("\n\tNumber of iterations\t\t\t %d", NIterations);
	printf("\n\tPercentage of Improvement\t\t\t %f", ZMin);
	printf("\n\tTeta\t\t\t %f", teta);
	printf("\n\tACO time limit\t\t %d", timeLimitACO);
	printf("\n\tRMP time limit\t\t %d", timeLimitRMP);
	printf("\n\tMIP time limit\t\t %d", timeLimitMIP);

}

void obtemValoresParametros() {

	int opcao = 14;

	while (opcao > 0) {

		printf("\nType:\n");
		printf("\t0 - To return.\n");
		printf("\t1 - To change QAnt.\n");
		printf("\t2 - To change Alpha.\n");
		printf("\t3 - To change Beta.\n");
		printf("\t4 - To change evaporation proportion.\n");
		printf("\t5 - To change addition proportion.\n");
		printf("\t6 - To change pheromone value.\n");
		printf("\t7 - To change QNZao.\n");
		printf("\t8 - To change Number of Iterations.\n");
		printf("\t9 - To change Percentage of Improvement.\n");
		printf("\t10 - To change Teta.\n");
		printf("\t11 - To change ACO time limit.\n");
		printf("\t12 - To change RMP Time Limit.\n");
		printf("\t13 - To change Total Limit.\n");
		scanf("%d", &opcao);

		system("cls");

		switch (opcao) {

		case 0:
			break;
		case 1:
			printf("\n\tType QAnt (Ex: 0.3)\n\t");
			scanf("%f", &QAnt);
			break;

		case 2:
			printf("\n\tType Alpha value (Ex: 1.0)\n\t");
			scanf("%f", &alfa);
			break;

		case 3:
			printf("\n\tType Beta value (Ex: 1.0)\n\t");
			scanf("%f", &beta);
			break;

		case 4:
			printf("\n\tType Evaporation proportion (Ex: 0.05)\n\t");
			scanf("%f", &evaporation);
			break;

		case 5:
			printf("\n\tType addition proportion (Ex: 0.05)\n\t");
			scanf("%f", &QIncrease);
			break;

		case 6:
			printf("\n\tType pheromone value (Ex: 1.0)\n\t");
			scanf("%f", &pheromone);
			break;

		case 7:
			printf("\n\tType QNZao value (Ex: 0.5)\n\t");
			scanf("%f", &QNZao);
			break;

		case 8:
			printf("\n\tType Number of Iterations (Ex: 20)\n\t");
			scanf("%d", &NIterations);
			break;

		case 9:
			printf("\n\tType Percentage of Improvement (Ex: 0.2)\n\t");
			scanf("%f", &ZMin);
			break;

		case 10:
			printf("\n\tType Teta (Ex: 0.75)\n\t");
			scanf("%f", &teta);
			break;

		case 11:
			printf("\n\tType ACO time limit in seconds(Ex: 10)\n\t");
			scanf("%d", &timeLimitACO);
			break;

		case 12:
			printf("\n\tType RMP time limit in seconds(Ex: 10)\n\t");
			scanf("%d", &timeLimitRMP);
			break;

		case 13:
			printf("\n\tType RMP (CG) time limit in seconds (Ex: 10)\n\t");
			scanf("%d", &timeLimitMIP);
			break;
		default:
			printf("\nPlease type a valid option\n");
			break;

		}

		if (opcao > 0 && opcao < 14) printf("\nAfter all changes, please check if they are processed correctly!\n");

	}

}


void importaParametros() {


	FILE* file;
	file = fopen("parametros-BPP-0-0-3.txt", "r");

	fscanf(file, "%f", &QAnt);
	fscanf(file, "%f", &alfa);
	fscanf(file, "%f", &beta);
	fscanf(file, "%f", &evaporation);
	fscanf(file, "%f", &QIncrease);
	fscanf(file, "%f", &pheromone);
	fscanf(file, "%f", &QNZao);
	fscanf(file, "%d", &NIterations);
	fscanf(file, "%f", &ZMin);
	fscanf(file, "%f", &teta);
	fscanf(file, "%d", &timeLimitACO);
	fscanf(file, "%d", &timeLimitRMP);
	fscanf(file, "%d", &timeLimitMIP);

	fclose(file);

}


void updateParams() {

	FILE* pont_arqu;

	pont_arqu = fopen("parametros-BPP-0-0-3.txt", "w");

	fprintf(pont_arqu, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%d\t%d\t%d\n", QAnt, alfa, beta, evaporation, QIncrease, pheromone, QNZao, NIterations, ZMin, teta, timeLimitACO, timeLimitRMP, timeLimitMIP);

	fclose(pont_arqu);

	system("cls");

	printParams();

}


//###################################### procedimentos MIP #####################################################

void alocaMemoriaRestricao(double** rhs, char** sense, char** senseMIP, int n) {

	*rhs = (double*)malloc(n * sizeof(double)); //reserva memória para o lado direito da restrição (double)

	//*sense = (char*)malloc(n * sizeof(char)); //reserva memória para o sinal das restrições (char)

	*senseMIP = (char*)malloc(n * sizeof(char)); //reserva memória para o sinal das restrições (char)

}

void criaRestricoesClass9(double* rhs, char* sense, char* senseMIP, int n, int* qtdTipos) {

	int i;

	for (i = 0; i < n; i++) {

		rhs[i] = qtdTipos[i]; //armazena o lado direito da restrição

		//sense[i] = GRB_GREATER_EQUAL; //armazena o sinal da restri��o

		senseMIP[i] = GRB_GREATER_EQUAL;

	}

}

void criaRestricoes(double* rhs, char* sense, char* senseMIP, int n) {

	int i;

	for (i = 0; i < n; i++) {

		rhs[i] = 1.0; //armazena o lado direito da restrição

		//sense[i] = GRB_GREATER_EQUAL; //armazena o sinal da restri��o

		senseMIP[i] = GRB_GREATER_EQUAL;

	}

}

void alocaMemoriaVariaveis(double** cMIPStart, int** cbeg, int** clen, int** cind, double** cval, double** lb, double** obj, char** ctype, char** ctypeMIP, int m, int nz) {

	*cbeg = (int*)malloc(m * sizeof(int)); //armazena a quantidade de variáveis (inteiro)

	*clen = (int*)malloc(m * sizeof(int)); //armazena a quantidade variáveis na restrição (inteiro)

	*cind = (int*)malloc(nz * sizeof(int)); //armazena os índices da restrição (inteiro)

	*cval = (double*)malloc(nz * sizeof(double)); //armazena os valores dos multiplicadores das variáveis (double)

	*lb = (double*)malloc(m * sizeof(double)); //armazena o limite inferior das variáveis (double)

	*obj = (double*)malloc(m * sizeof(double)); //armazena o multiplicador das variáveis na função objetivo (double)

	//*ctype = (char*)malloc(m * sizeof(char)); //armazena o tipo de variaval (Integer, Binary, Continuous)

	*ctypeMIP = (char*)malloc(m * sizeof(char)); //armazena o tipo de variaval (Integer, Binary, Continuous)

	*cMIPStart = (double*)calloc(m, sizeof(double));

}

void preencheVariaveisClass9(double* cMIPStart, NoBin** ultimaColuna, NoSoluction* soluctions, int* cbeg, int* cind, double* cval, int* clen, char* ctype, char* ctypeMIP, int quantidade, int** itens, int totalTipo) {

	int* usedItens = malloc(totalTipo * sizeof(int));

	int* indicesBestStart = NULL, bestActual = 0, totalIndiceBest = 0;

	int nz = 0, initialNz, actualNz, totalNZ, tipo;

	int i = 0, k = 0, j = 0;

	int* indiceBestStart = NULL;

	float valorMIPStart = soluctions->soluction.utilization;

	NoSoluction* auxSoluction = soluctions;

	NoBin* auxBin, * binAnterior = NULL;

	NoPack* auxPack;

	while (auxSoluction) {

		bestActual = 0;


		if (indicesBestStart == NULL || auxSoluction->soluction.utilization < valorMIPStart) {

			if (indicesBestStart) free(indicesBestStart);

			totalIndiceBest = auxSoluction->soluction.value;

			indicesBestStart = calloc(auxSoluction->soluction.value, sizeof(int));

			bestActual = 1;

			valorMIPStart = auxSoluction->soluction.utilization;

			j = 0;

		}

		auxBin = auxSoluction->soluction.bins;

		while (auxBin) {

			cMIPStart[i] = 0.0;

			initialNz = nz;

			auxBin->bin.idtColuna = i;

			cbeg[i] = nz;

			auxPack = auxBin->bin.conteudo;

			for (k = 0; k < totalTipo; k++) {

				usedItens[k] = -1;
			}

			totalNZ = 0;

			while (auxPack) {

				tipo = itens[auxPack->pack.id][5];

				if (usedItens[tipo] == -1) {

					usedItens[tipo] = nz + totalNZ;

					cval[nz + totalNZ] = 0.0;

					cind[nz + totalNZ] = tipo;

					totalNZ++;

				}

				actualNz = usedItens[tipo];

				cval[actualNz] = cval[actualNz] + 1.0;

				auxPack = auxPack->proximo;

			}

			nz = nz + totalNZ;

			if (bestActual == 1) {//marca a variavel no MIP Start

				indicesBestStart[j] = i;

				j++;

			}

			//ctype[i] = GRB_CONTINUOUS;

			ctypeMIP[i] = GRB_INTEGER;

			clen[i] = nz - cbeg[i];

			i++;

			if (binAnterior) {

				binAnterior->proximoColuna = auxBin;

			}

			binAnterior = auxBin;

			auxBin = auxBin->proximo;

		}

		auxSoluction = auxSoluction->proximo;

	}

	for (k = 0; k < totalIndiceBest; k++) {

		cMIPStart[indicesBestStart[k]] = (double)1.0;

	}

	free(indicesBestStart);

	free(usedItens);

	*ultimaColuna = binAnterior;

}


void preencheVariaveis(double* cMIPStart, NoBin** ultimaColuna, NoSoluction* soluctions, int* cbeg, int* cind, double* cval, int* clen, char* ctype, char* ctypeMIP, int quantidade) {

	int nz = 0, i = 0;

	int* indicesBestStart = NULL, bestActual = 0, totalIndiceBest = 0;

	int k = 0, j = 0;

	int* indiceBestStart = NULL;

	float valorMIPStart = soluctions->soluction.utilization;

	NoSoluction* auxSoluction = soluctions;

	NoBin* auxBin, * binAnterior = NULL;

	NoPack* auxPack;

	while (auxSoluction) {

		bestActual = 0;

		if (indicesBestStart == NULL || auxSoluction->soluction.utilization < valorMIPStart) {

			if (indicesBestStart) free(indicesBestStart);

			totalIndiceBest = auxSoluction->soluction.value;

			indicesBestStart = calloc(auxSoluction->soluction.value, sizeof(int));

			valorMIPStart = auxSoluction->soluction.utilization;

			bestActual = 1;

		}

		auxBin = auxSoluction->soluction.bins;

		j = 0;

		while (auxBin) {

			auxBin->bin.idtColuna = i;

			cbeg[i] = nz;

			auxPack = auxBin->bin.conteudo;

			while (auxPack) {

				cind[nz] = auxPack->pack.id;

				cval[nz] = 1.0;

				nz++;

				auxPack = auxPack->proximo;

			}

			if (bestActual == 1) {//marca a variavel no MIP Start

				indicesBestStart[j] = i;

				j++;

			}

			ctypeMIP[i] = GRB_BINARY;

			//ctype[i] = GRB_CONTINUOUS;

			clen[i] = nz - cbeg[i];

			i++;

			if (binAnterior) {

				binAnterior->proximoColuna = auxBin;

			}

			binAnterior = auxBin;

			auxBin = auxBin->proximo;

		}

		auxSoluction = auxSoluction->proximo;

	}

	for (k = 0; k < totalIndiceBest; k++) {

		cMIPStart[indicesBestStart[k]] = (double)1.0;

	}

	free(indicesBestStart);

	*ultimaColuna = binAnterior;

}

void criaFuncaoObjetivo(double* lb, double* obj, int m) {

	int j;

	for (j = 0; j < m; j++) {
		//for x variable

		lb[j] = 0.0; //lower bound da variável x

		obj[j] = 1.0; //valor do objetivo = c

	}

}

int retornaQuantidadeNaoZerosBin(NoPack* conteudo, int totalTipo, int** itens) {

	NoPack* auxPack = conteudo;

	int nz = 0;

	int* usedItens = calloc(totalTipo, sizeof(int));

	while (auxPack) {

		if (usedItens[itens[auxPack->pack.id][5]] == 0) {

			usedItens[itens[auxPack->pack.id][5]] = 1;

			nz++;

		}

		auxPack = auxPack->proximo;

	}

	free(usedItens);

	return nz;

}

int retornaQuantidadeNaoZerosSoluction(NoSoluction* soluction, int totalTipo, int** itens) {

	NoBin* auxBin = soluction->soluction.bins;

	int nz = 0;

	while (auxBin) {


		nz = nz + retornaQuantidadeNaoZerosBin(auxBin->bin.conteudo, totalTipo, itens);

		auxBin = auxBin->proximo;

	}

	return nz;
}

int montaModeloClass9(GRBenv* env, GRBmodel** model, GRBmodel** modelMIP, NoBin** ultimaColuna, NoSoluction* geralSoluction, int* m, int totalTipo, int* qtdTipos, int** itens) {

	int error = 0;
	int nz = 0;
	int n = totalTipo;
	int* cbeg = NULL; //indica onde a restrição se inicia no array cind
	int* clen = NULL;//indica quantos índices não nulos existem na restrição
	int* cind = NULL;//contém os indices da variável na restrição
	double* cval = NULL;//indica o multiplicador da variavel na restrição
	double* rhs = NULL; //lado direito da restrição
	char* sense = NULL; // sinal da restrição
	char* senseMIP = NULL;
	double* lb = NULL; //lower bound das variáveis
	double* obj = NULL; //armazena os valores da função objetivo
	char* ctype = NULL; //tipo da vari�vel
	char* ctypeMIP = NULL; //tipo da variável
	double* cMIPStart = NULL;

	NoSoluction* auxSoluction = geralSoluction;
	NoBin* topoBin = NULL;

	alocaMemoriaRestricao(&rhs, &sense, &senseMIP, n);


	criaRestricoesClass9(rhs, sense, senseMIP, n, qtdTipos);

	while (auxSoluction) {

		*m = *m + auxSoluction->soluction.value; // quantidade de variáveis

		nz = nz + retornaQuantidadeNaoZerosSoluction(auxSoluction, totalTipo, itens); //quantidade de indices nao zero nas restricoes

		auxSoluction = auxSoluction->proximo;

	}

	alocaMemoriaVariaveis(&cMIPStart, &cbeg, &clen, &cind, &cval, &lb, &obj, &ctype, &ctypeMIP, *m, nz);

	criaFuncaoObjetivo(lb, obj, *m);

	preencheVariaveisClass9(cMIPStart, ultimaColuna, geralSoluction, cbeg, cind, cval, clen, ctype, ctypeMIP, *m, itens, totalTipo);

	/*error = GRBloadmodel(env, model, "RMP", *m, n,
		GRB_MINIMIZE, 0.0, obj, sense, rhs,
		cbeg, clen, cind, cval, lb, NULL,
		ctype, NULL, NULL);
	if (error) {printf("\nErro ao montar o modelo RMP (2)\n");  return error;}
	*/
	error = GRBloadmodel(env, modelMIP, "MIP", *m, n,
		GRB_MINIMIZE, 0.0, obj, senseMIP, rhs,
		cbeg, clen, cind, cval, lb, NULL,
		ctypeMIP, NULL, NULL);
	if (error) { printf("\nErro ao montar o modelo MIP (2)\n");   return error; }


	error = GRBsetdblattrarray(*modelMIP, "Start", 0, *m, cMIPStart);
	if (error) { printf("\nErro ao cadastrar o MIP Start (2)\n");   return error; }

	return error;

}

int montaModelo(GRBenv* env, GRBmodel** model, GRBmodel** modelMIP, NoBin** ultimaColuna, NoSoluction* geralSoluction, int quantidade, int* m) {

	int error = 0;
	int nz = 0;
	int n = quantidade;
	int* cbeg = NULL; //indica onde a restrição se inicia no array cind
	int* clen = NULL;//indica quantos índices não nulos existem na restrição
	int* cind = NULL;//contém os indices da variável na restrição
	double* cval = NULL;//indica o multiplicador da variavel na restrição
	double* rhs = NULL; //lado direito da restrição
	char* sense = NULL; // sinal da restrição
	char* senseMIP = NULL;
	double* lb = NULL; //lower bound das variáveis
	double* obj = NULL; //armazena os valores da função objetivo
	char* ctype = NULL; //tipo da vari�vel
	char* ctypeMIP = NULL; //tipo da variável
	double* cMIPStart = NULL;

	NoSoluction* auxSoluction = geralSoluction;
	NoBin* topoBin = NULL;

	alocaMemoriaRestricao(&rhs, &sense, &senseMIP, n);

	criaRestricoes(rhs, sense, senseMIP, n);

	while (auxSoluction) {

		*m = *m + auxSoluction->soluction.value; // quantidade de variáveis

		nz = nz + n; //quantidade de indices nao zero nas restricoes

		auxSoluction = auxSoluction->proximo;

	}

	alocaMemoriaVariaveis(&cMIPStart, &cbeg, &clen, &cind, &cval, &lb, &obj, &ctype, &ctypeMIP, *m, nz);

	criaFuncaoObjetivo(lb, obj, *m);

	preencheVariaveis(cMIPStart, ultimaColuna, geralSoluction, cbeg, cind, cval, clen, ctype, ctypeMIP, *m);

	/*error = GRBloadmodel(env, model, "RMP", *m, n,
		GRB_MINIMIZE, 0.0, obj, sense, rhs,
		cbeg, clen, cind, cval, lb, NULL,
		ctype, NULL, NULL);
	if (error) { printf("\nErro ao montar o modelo RMP (1)\n");  return error; }
	*/
	error = GRBloadmodel(env, modelMIP, "MIP", *m, n,
		GRB_MINIMIZE, 0.0, obj, senseMIP, rhs,
		cbeg, clen, cind, cval, lb, NULL,
		ctypeMIP, NULL, NULL);
	if (error) { printf("\nErro ao montar o modelo MIP (1)\n");   return error; }

	error = GRBsetdblattrarray(*modelMIP, "Start", 0, *m, cMIPStart);
	if (error) { printf("\nErro ao cadastrar o MIP Start (1)\n");   return error; }

	return error;

}

int obtemResultado(GRBmodel* model, double** valorVariaveis, int m) {

	int erro;

	*valorVariaveis = (double*)malloc(m * sizeof(double));

	erro = GRBgetdblattrarray(model, "X", 0, m, *valorVariaveis);
	if (erro) printf("\nErro ao obter o valor das variaveis\n");

	return erro;
}

void salvaSolucao(NoSoluction** columnSoluction, NoBin* primeiraColuna, double* valorVariaveis, int m, int valor, int totalTime) {

	NoBin* auxBin = primeiraColuna, * topoBin = NULL, * deleteBin = NULL;


	int j;

	for (j = 0; j < m; j++) {

		if (valorVariaveis[j] > 0.0) {

			auxBin->bin.idt = (int)valorVariaveis[j];
			auxBin->proximo = topoBin;
			topoBin = auxBin;
			deleteBin = NULL;

		}
		else {

			deleteBin = auxBin;
		}


		auxBin = auxBin->proximoColuna;

		if (deleteBin) {

			deleteBin->proximo = NULL;
			freeMemoryBin(deleteBin);

		}




	}

	*columnSoluction = empilharSoluction(*columnSoluction, topoBin, 0, valor, 0.0, totalTime);

}

void excluiItensDuplicados(NoSoluction* columnSoluction, int** itens, int quantidade) {

	int* usedItens = (int*)calloc(quantidade, sizeof(int));

	NoBin* auxBin = columnSoluction->soluction.bins;

	NoPack* auxPack = NULL;

	while (auxBin) {

		auxPack = auxBin->bin.conteudo;

		while (auxPack) {

			if (usedItens[auxPack->pack.id] == 0) {// item ainda não utilizado

				usedItens[auxPack->pack.id] = 1;

			}
			else {

				deletarPack(auxBin, auxPack->pack.id, itens[auxPack->pack.id][1] * itens[auxPack->pack.id][2] * itens[auxPack->pack.id][3]);

			}

			auxPack = auxPack->proximo;

		}

		auxBin = auxBin->proximo;
	}


	free(usedItens);

}

int selectBestCombination(int** combinacoes, float* probabilityCombination, int quantidadeCombinacoes) {

	int k, bestCombination = 0;
	float bestCombinationValue = 0.0;

	for (k = 0; k < quantidadeCombinacoes; k++) {

		if (bestCombinationValue < probabilityCombination[k]) {

			bestCombinationValue = probabilityCombination[k];

			bestCombination = k;

		}

	}

	o = combinacoes[bestCombination][1];
	sr = combinacoes[bestCombination][2];

	return bestCombination;

}

void  valoresIniciaisDeOrdemCombinado(int** itens, double* PiValue, float* nextStep, int* ordem, int quantidade, int classe, float* heuristica) {

	int j;

	for (j = 0; j < quantidade; j++) {

		if (classe == 9) {

			ordem[j] = j;
			nextStep[j] = (float)PiValue[itens[j][5]] * heuristica[j];

		}
		else {

			ordem[j] = j;
			nextStep[j] = (float)PiValue[j] * heuristica[j];

		}


	}

}


void  valoresIniciaisDeOrdem(int** itens, double* PiValue, float* nextStep, int* ordem, int quantidade, int classe) {

	int j;

	for (j = 0; j < quantidade; j++) {

		if (classe == 9) {

			ordem[j] = j;
			nextStep[j] = (float)PiValue[itens[j][5]];

		}
		else {

			ordem[j] = j;
			nextStep[j] = (float)PiValue[j];

		}


	}

}

/*
void testeOrdenamento(int* ordem, float* nextStep, int quantidade) {

	int i;
	float anterior = nextStep[ordem[0]];

	for (i = 0; i < quantidade; i++) {

		if (nextStep[ordem[i]] > anterior) {

			printf("\nErro no ordenamento! Atual %f, maior que anterior %f\n");

			system("pause");

		}

		anterior = nextStep[ordem[i]];

	}

	printf("\nPassou no teste Ordenamento\n");

	system("pause");

}*/


NoBin* CriaNovosPadroes(int combinacao, int** itens, float** heuristicas, int* ordem, double* PiValue, int quantidade, int tamanho[], int classe, float* nextStep) {

	NoBin* auxBin = NULL, * auxBin1 = NULL;
	NoBin* topoBin = NULL;
	int k, l = 0;

	valoresIniciaisDeOrdem(itens, PiValue, nextStep, ordem, quantidade, classe);

	//printf("\nAtribui valores iniciais de ordem\n");
	//system("pause");

	quicksortDecrease(0, quantidade - 1, nextStep, ordem);

	//printf("\nOrdenou itens\n");
	//system("pause");

//	for (k = 0; k < 4; k++) {

		//printf("\nDefinido o\n");
		//system("pause");

//		o = combinacoes[k][1];

		//printf("\nDefinido o\n");
		//system("pause");

//		sr = combinacoes[k][2];

		//printf("\nDefinido sr\n");
		//system("pause");

	auxBin = binPack(itens, ordem, quantidade, tamanho, PiValue, classe);

	//printf("\nSolucao com combinacao %i\n", k);
	//system("pause");

//		auxBin1 = auxBin;

//		while (auxBin1->proximo) {

//			auxBin1 = auxBin1->proximo;

//		}

//		auxBin1->proximo = topoBin;
	topoBin = auxBin;

	//	}

		//printf("\nCriou padrao combinado apenas com dual\n");
		//system("pause");

	//	while (l < 10) {

	valoresIniciaisDeOrdemCombinado(itens, PiValue, nextStep, ordem, quantidade, classe, heuristicas[combinacao]);

	quicksortDecrease(0, quantidade - 1, nextStep, ordem);

	//		for (k = 0; k < 4; k++) {

	//			o = combinacoes[k][1];

	//			sr = combinacoes[k][2];

	auxBin = binPack(itens, ordem, quantidade, tamanho, PiValue, classe);

	auxBin1 = auxBin;

	while (auxBin1->proximo) {

		auxBin1 = auxBin1->proximo;

	}

	auxBin1->proximo = topoBin;
	topoBin = auxBin;

	//		}

	//		l = l + 2;
	//	}


	//	printf("\nCriou padrao combinado dual e heuristica\n");
	//	system("pause");

	return topoBin;

}

void alocaMemoriaUmaVariavel(int** cind, double** cval, int quantidade) {

	*cind = (int*)malloc(quantidade * sizeof(int)); //armazena os �ndices da restri��o (inteiro)

	*cval = (double*)malloc(quantidade * sizeof(double)); //armazena os valores dos multiplicadores das vari�veis (double)

}

void zeraQtdTiposBin(int* qtdTiposBin, int totalTipos) {

	int j;

	for (j = 0; j < totalTipos; j++) {

		qtdTiposBin[j] = 0;

	}

}

int preencheUmaVariavel(int* cind, double* cval, NoBin* auxBin) {

	int nz = 0;

	NoPack* auxPack = auxBin->bin.conteudo;

	while (auxPack) {

		cind[nz] = auxPack->pack.id;

		cval[nz] = 1.0;

		nz++;

		auxPack = auxPack->proximo;
	}

	return nz;

}

int preencheUmaVariavelClass9(int* cind, double* cval, int qtdTipos, int* totalTiposNoBin) {

	int nz = 0, k;

	for (k = 0; k < qtdTipos; k++) {

		if (totalTiposNoBin[k] > 0) {

			cind[nz] = k;

			cval[nz] = totalTiposNoBin[k];

			nz++;
		}

	}

	return nz;
}


NoBin* insereColuna(NoBin* colunaAInserir, NoBin* ultimaColunaInserida, int idColuna) {

	colunaAInserir->bin.idtColuna = idColuna;
	colunaAInserir->proximoColuna = NULL;
	ultimaColunaInserida->proximoColuna = colunaAInserir;
	ultimaColunaInserida = colunaAInserir;

	return ultimaColunaInserida;


}



int fixaUmaVariavel(GRBmodel* model, double* valorVariaveis, int m, int* custo_reduzido_negativo) {

	int k, id = -1, error = 0, arredondado = 0;
	double maiorValor = 0.0, fractionalPart, integerPart;

	for (k = 0; k < m; k++) {

		if (valorVariaveis[k] > 0.0 && floor(valorVariaveis[k]) < ceil(valorVariaveis[k])) {

			fractionalPart = modf(valorVariaveis[k], &integerPart);

			if (fractionalPart >= teta) {

				error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, k, integerPart + 1.0);
				if (error) { printf("\nErro ao atribuir Lower Bound a variavel(1)\n"); return error; }

				if (arredondado == 0) {
					*custo_reduzido_negativo = 1;
					arredondado = 1;
				}

			}
			else if (arredondado == 0 && maiorValor < fractionalPart) {

				id = k;
				maiorValor = fractionalPart;

			}

		}

	}

	if (id > -1 && arredondado == 0) {


		fractionalPart = modf(valorVariaveis[id], &integerPart);

		*custo_reduzido_negativo = 1;
		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, id, integerPart + 1.0);

		//printf("\n\nArredondada variavel %lf para %lf\n", valorVariaveis[id], integerPart + 1.0);

	}

	if (error) printf("\nErro ao atribuir Lower Bound a variavel(2)\n");

	return error;

}


int verificaSeEhInteiro(GRBmodel* model, int m) {

	int i;

	double* values = (double*)malloc(m * sizeof(double));

	GRBgetdblattrarray(model, "X", 0, m, values);

	for (i = 0; i < m; i++) {

		if (values[i] != 0.0 && values[i] != 1.0) {

			return 0;

		}

	}

	free(values);

	return 1;

}

//##################################### testes RMP ################################################
/*
void testeIntegridadeColunas(NoBin* primeiraColuna, NoBin* ultimaColuna, int m) {

	int i=0;
	NoBin* auxBin = primeiraColuna, *binAnterior = primeiraColuna;

	while (auxBin) {

		i++;

		if (i > m) {

			printf("\nErro! Variavel m: %d, bin atual: %d\n", m, i);
			system("pause");

		}

		binAnterior = auxBin;

		auxBin = auxBin->proximoColuna;
	}

	if(i<m) {

		printf("\nErro! Variavel m: %d, bin atual: %d\n", m, i);
		system("pause");

	}

	if(ultimaColuna != binAnterior) {

		printf("\nErro! Ultima Coluna: %p, ultima Bin: %p\n", ultimaColuna, binAnterior);
		system("pause");

	}

	printf("\nPassou no teste Integridade Colunas\n");

	system("pause");

}

void testeQtdTiposBin(int* qtdTiposBin, int totalTipos, int totalTiposBin, NoPack* topoPack, int **itens) {

	NoPack* auxPack = topoPack;

	int* referencia = calloc(totalTipos, sizeof(int));

	int nz = 0;

	while (auxPack) {

		if (referencia[itens[auxPack->pack.id][5]] == 0) nz++;

		referencia[itens[auxPack->pack.id][5]] = referencia[itens[auxPack->pack.id][5]] + 1;

		if (referencia[itens[auxPack->pack.id][5]] > qtdTiposBin[itens[auxPack->pack.id][5]]) {

			printf("\nErro! Referencia: %d, qtdTiposBin: %d\n", referencia[itens[auxPack->pack.id][5]], qtdTiposBin[itens[auxPack->pack.id][5]]);
			system("pause");

		}

		auxPack = auxPack->proximo;

	}

	if (nz != totalTiposBin) {

		printf("\nErro! nz: %d, totalTiposBin: %d\n", nz, totalTiposBin);
		system("pause");

	}

	for (nz = 0; nz < totalTipos; nz++) {

		if (referencia[nz] != qtdTiposBin[nz]) {

			printf("\nErro! tipo %: referencia %d, qtdTiposBin: %d\n", nz, referencia[nz], qtdTiposBin[nz]);
			system("pause");

		}

	}


	free(referencia);

	printf("\nPassou no teste Quantidade Tipos na Bin\n");

	system("pause");

}
*/

void printLog(int passo) {

	FILE* pont_arqu;

	pont_arqu = fopen("log.txt", "a");

	fprintf(pont_arqu, " %d", passo);

	if (passo == 42) fprintf(pont_arqu, "\n");

	fclose(pont_arqu);

}

int defineSeed(int classe, int quantidade, int instancia) {

	int semente = 0, total, contador, classeAtual, quantidadeAtual, instanciaAtual, temp;

	char texto[30];
	char texto1[3];
	char texto2[] = ".txt";

	strcpy(texto, "seeds/default-seeds-");

	sprintf(texto1, "%d", classe);

	strcat(texto, texto1);
	strcat(texto, texto2);

	FILE* file;
	file = fopen(texto, "r");

	fscanf(file, "%i", &total);

	if (classe != 9) {

		for (contador = 0; contador < total; contador++) {

			fscanf(file, "%i", &quantidadeAtual);
			fscanf(file, "%i", &instanciaAtual);
			fscanf(file, "%i", &temp);

			if (quantidadeAtual == quantidade && instanciaAtual == instancia) {

				semente = temp;
				contador = total;

			}

		}

	}
	else {

		for (contador = 0; contador < total; contador++) {

			fscanf(file, "%i", &instanciaAtual);
			fscanf(file, "%i", &temp);

			if (instanciaAtual == instancia) {

				semente = temp;
				contador = total;

			}

		}
	}

	fclose(file);

	if (semente == 0) {

		printf("Erro!Semente nao encontrada");
		system("pause");
	}

	return semente;
}


void defineSettings(int classe, int quantidade) {

	int actualACOTimeLimit = 0, actualRMPTimeLimit = 0, actualMIPTimeLimit = 0, total, contador, classeAtual, quantidadeAtual, temp;

	char texto[30];

	FILE* file;

	if (classe < 9) {

		strcpy(texto, "settings/times-1-8.txt");

		file = fopen(texto, "r");

		fscanf(file, "%i", &total);

		for (contador = 0; contador < total; contador++) {

			fscanf(file, "%i", &classeAtual);
			fscanf(file, "%i", &quantidadeAtual);
			fscanf(file, "%i", &temp);
			fscanf(file, "%i", &actualRMPTimeLimit);
			fscanf(file, "%i", &actualMIPTimeLimit);

			if (quantidadeAtual == quantidade && classeAtual == classe) {

				actualACOTimeLimit = temp;
				contador = total;

			}

		}

	}
	else {

		strcpy(texto, "settings/times-9-10.txt");

		file = fopen(texto, "r");

		fscanf(file, "%i", &total);

		for (contador = 0; contador < total; contador++) {

			fscanf(file, "%i", &quantidadeAtual);
			fscanf(file, "%i", &temp);
			fscanf(file, "%i", &actualRMPTimeLimit);
			fscanf(file, "%i", &actualMIPTimeLimit);

			if (quantidadeAtual == classe || quantidadeAtual == quantidade) {

				actualACOTimeLimit = temp;
				contador = total;

			}

		}

	}

	fclose(file);

	if (actualACOTimeLimit == 0) {

		printf("Erro!Nao encontrado os tempos limites em settings");
		system("pause");

	}
	else {

		timeLimitACO = actualACOTimeLimit;
		timeLimitRMP = actualRMPTimeLimit;
		timeLimitMIP = actualMIPTimeLimit;

	}

}

void defineQAnt(int quantidade) {

	int total, contador, quantidadeAtual;
	float temp, actualQAnt = 0.0;

	char texto[30];

	FILE* file;

	strcpy(texto, "settings/QAnt-10.txt");

	file = fopen(texto, "r");

	fscanf(file, "%i", &total);

	for (contador = 0; contador < total; contador++) {

		fscanf(file, "%i", &quantidadeAtual);
		fscanf(file, "%f", &temp);

		if (quantidadeAtual == quantidade) {

			actualQAnt = (float)temp;
			contador = total;

		}

	}

	fclose(file);

	if (actualQAnt == 0.0) {

		printf("Erro!Nao encontrado os tempos limites em settings");
		system("pause");

	}
	else {

		QAnt = actualQAnt;

	}

}

int main() {

	//#####################variaveis gerais#############################################

	int opcao, tamanho[3], classe = 0, quantidade, instancia, totalTime, initialTime, tempoUltimaMelhoria, espera;
	unsigned int semente = (unsigned)time(NULL);
	int** itens;
	int k;
	int continuousLowerBound;

	int ants, actualAnt, bestCombination;
	int quantidadeCombinacoes;
	int quantidadeOrdens;
	float utilization, decrease, totalProbabilityCombination;
	int* ordem, * bestOrder;
	float** heuristicas;
	float** pheromones;
	int **combinacoes;
	float *denominator, *probabilityCombination;

	//#####################variaveis Solucao#####################

	int idSoluction, erro;
	NoBin* topoBin, * primeiraColuna = NULL, * ultimaColuna = NULL, * auxBin = NULL;
	NoPack* auxPack = NULL;
	NoSoluction* geralSoluction, * bestOfBest, * bestIteration, * coveringSoluction;
	NoSoluction* columnSoluction = NULL, * roudingSoluction = NULL;

	//#######################variaveis RMP######################
	GRBenv* env = NULL; //cria��o do ambiente
	GRBmodel* model = NULL;
	int idColuna = 0; //id da �ltima coluna adicionada
	int custo_reduzido_negativo = 1; // se h� custo reduzido negativo
	int integer = 0;  //se a solu��o � inteira
	int totalTipos;
	double reducedCost; //valor do custo reduzido do padr�o atual
	double* PiValue;
	int actualIteration = 0;
	double referenceImprovement;
	int totalTiposBin, * qtdTipos = NULL, * qtdTiposBin = NULL;
	float* nextStep;

	//#######################variaveis MIP######################
	GRBmodel* modelMIP = NULL; //cria��o do modelo
	double* valorVariaveis = NULL;
	int  nz = 0, i, j; //nz n�o zeros, i = x, j = contador
	int status;
	double objval;
	int parcialTime;
	int error = 0, m = 0, n, nova_m = 0;

	error = GRBloadenv(&env, NULL);
	if (error) { printf("\nErro ao montar o ambiente\n"); goto QUIT; }

	//error = GRBgetdblparam(env, GRB_DBL_PAR_INTFEASTOL, &gurobiTolerance);
	//if (error) { printf("\nErro ao obter o valor da IntFeasTol\n"); goto QUIT; }

	importaParametros();

	printf("\tBPP Version 0.1.6\t\tMethods EP\n\tDetails: Insert News Combinations");
	printf("\n\tAuthor: Daniel Bento Maia\tDate: 11/11/2025\n");

	if (classe == 9) {

		quantidadeCombinacoes = 10;
		quantidadeOrdens = 5;

	}
	else {

		quantidadeCombinacoes = 7;
		quantidadeOrdens = 4;

	}

	/*	printf("\n\n\tTo change the parameters type 2; to solve one problem, type 1\n\t");
		scanf("%d", &opcao);

		if (opcao == 2) {

			obtemValoresParametros();

			updateParams();

			opcao = 1;
		}
	*/

	printf("\n\n\tPlease type class, quantity of itens and instance:");
	printf("\n\tFor example, to solve file 1_50_3.txt type 1 50 3\n\n\t");
	scanf("%d %d %d", &classe, &quantidade, &instancia);


	denominator = malloc(quantidadeCombinacoes * sizeof(float));
	probabilityCombination = malloc(quantidadeCombinacoes * sizeof(float));

	combinacoes = malloc(quantidadeCombinacoes * sizeof(int*));
	pheromones = malloc(quantidadeCombinacoes * sizeof(float*));
	heuristicas = malloc(quantidadeOrdens * sizeof(float*));

	defineSettings(classe, quantidade);

	espera = (timeLimitACO * 0.15) + 1;

	if (classe == 10) defineQAnt(quantidade);

	itens = malloc(quantidade * sizeof(int*));
	ordem = (int*)malloc(quantidade * sizeof(int));
	bestOrder = (int*)malloc(quantidade * sizeof(int));
	nextStep = malloc(quantidade * sizeof(float));

	for (k = 0; k < quantidade; k++) {

		itens[k] = malloc(6 * sizeof(int));

	}

	for (k = 0; k < quantidadeOrdens; k++) {

		heuristicas[k] = malloc(quantidade * sizeof(float));

	}

	for (k = 0; k < quantidadeCombinacoes; k++) {

		pheromones[k] = malloc(quantidade * sizeof(float));
		combinacoes[k] = malloc(3 * sizeof(int));

	}

	preencheCombinacoes(combinacoes, quantidadeCombinacoes);

	semente = defineSeed(classe, quantidade, instancia);

	srand(semente);

	//	printParams();

	tamanhoInstancia(classe, tamanho);

	totalTime = 0;
	idSoluction = 0;
	topoBin = NULL;
	geralSoluction = NULL;
	bestIteration = NULL;
	bestOfBest = NULL;
	roudingSoluction = NULL;
	primeiraColuna = NULL;
	ultimaColuna = NULL;
	auxPack = NULL;
	columnSoluction = NULL;
	roudingSoluction = NULL;
	model = NULL;
	idColuna = 0;
	custo_reduzido_negativo = 1;
	integer = 0;
	actualIteration = 0;
	qtdTipos = NULL;
	qtdTiposBin = NULL;
	modelMIP = NULL;
	valorVariaveis = NULL;
	nz = 0;
	error = 0; m = 0; nova_m = 0;


	if (classe == 9) {

		totalTipos = importaArquivoClass9(itens, heuristicas, classe, quantidade, instancia, tamanho);

		PiValue = (double*)malloc(totalTipos * sizeof(double));

		qtdTipos = calloc(totalTipos, sizeof(int));

		for (k = 0; k < quantidade; k++) {

			qtdTipos[itens[k][5]] = qtdTipos[itens[k][5]] + 1;

		}

	}
	else {

		totalTipos = quantidade;

		importaArquivo(itens, heuristicas, classe, quantidade, instancia, tamanho);

		PiValue = (double*)malloc(quantidade * sizeof(double));

	}

	for (k = 0; k < totalTipos; k++) {

		PiValue[k] = 0.1;

	}

	printParams();

	initialTime = time(NULL);

	firstOrientation(itens, quantidade);

	NoProbability* topoProbability = NULL;

	topoProbability = montaPilhaProbabilidade(quantidade, pheromones, heuristicas, quantidadeCombinacoes, quantidadeOrdens);

	ants = ((float)quantidade) * QAnt;

	decrease = (float)pheromone / (float)quantidade;

	continuousLowerBound = volumeTotal / (tamanho[0] * tamanho[1] * tamanho[2]);

	totalProbabilityCombination = 0.0;

	for (k = 0; k < quantidadeCombinacoes; k++) {

		ordenaValores(heuristicas[combinacoes[k][0]], ordem, quantidade);

		o = combinacoes[k][1];

		sr = combinacoes[k][2];

		topoBin = binPack(itens, ordem, quantidade, tamanho, PiValue, classe);

		utilization = calcSmallerUtilization(topoBin, tamanho);

		geralSoluction = empilharSoluction(geralSoluction, topoBin, idSoluction, topoBin->bin.idt + 1, (float)topoBin->bin.idt + (float)1.0 + utilization, 0);

		idSoluction++;

		probabilityCombination[k] = ((float)continuousLowerBound) / geralSoluction->soluction.utilization;

		totalProbabilityCombination = totalProbabilityCombination + probabilityCombination[k];

		if (bestOfBest == NULL || bestOfBest->soluction.utilization > geralSoluction->soluction.utilization) {


			bestOfBest = geralSoluction;

		}

	}

	printf("\n\tSolucao %d\t\t%2f", bestOfBest->soluction.id, bestOfBest->soluction.utilization);

	calculaProbabilityCombination(probabilityCombination, totalProbabilityCombination, quantidadeCombinacoes);

	tempoUltimaMelhoria = time(NULL);

	//AQUI ACONTECE O ACO
	while (totalTime < timeLimitACO && time(NULL) - tempoUltimaMelhoria < espera) {

		actualAnt = 0;

		calcDenominator(combinacoes, heuristicas, pheromones, quantidade, topoProbability, denominator, quantidadeCombinacoes);

		bestIteration = NULL;

		while (actualAnt <= ants) {

			k = sorteiaProbabilidadeCombinacao(probabilityCombination, quantidadeCombinacoes);

			calcProbability(quantidade, denominator, topoProbability, ordem, k);

			topoBin = binPack(itens, ordem, quantidade, tamanho, PiValue, classe);

			utilization = calcSmallerUtilization(topoBin, tamanho);

			if (bestIteration == NULL || bestIteration->soluction.utilization > (float)topoBin->bin.idt + (float)1.0 + utilization) {

				if (bestIteration != NULL) {


					freeMemorySoluction(bestIteration);
					bestIteration = NULL;

				}

				bestIteration = empilharSoluction(bestIteration, topoBin, idSoluction, topoBin->bin.idt + 1, (float)topoBin->bin.idt + (float)1 + utilization, totalTime);

				idSoluction++;

				copyVector(ordem, bestOrder, quantidade);

				bestCombination = k;
			}
			else {

				freeMemoryBin(topoBin);

			}

			totalTime = time(NULL) - initialTime;

			actualAnt++;

		}

		bestIteration->proximo = geralSoluction;

		geralSoluction = bestIteration;

		totalProbabilityCombination = updateProbabilityCombination(probabilityCombination, totalProbabilityCombination, bestCombination, continuousLowerBound / bestIteration->soluction.utilization, quantidadeCombinacoes);

		updatePheromone(bestOrder, pheromones[k], quantidade, continuousLowerBound / bestIteration->soluction.utilization, decrease);

		if (bestIteration->soluction.utilization < bestOfBest->soluction.utilization) {

			bestOfBest = bestIteration;

			tempoUltimaMelhoria = time(NULL);

			printf("\n\tSolucao %d\t\t%2f", bestOfBest->soluction.id, bestOfBest->soluction.utilization);

		}


	}

	bestOfBest->soluction.time = time(NULL) - initialTime;

	primeiraColuna = geralSoluction->soluction.bins;

	bestCombination = selectBestCombination(combinacoes, probabilityCombination, quantidadeCombinacoes);

	referenceImprovement = (double)bestOfBest->soluction.utilization;

	//printLog(0);


	if (classe == 9) {
		error = montaModeloClass9(env, &model, &modelMIP, &ultimaColuna, geralSoluction, &m, totalTipos, qtdTipos, itens);
		if (error) goto QUIT;

		//printLog(1);

	}
	else {
		error = montaModelo(env, &model, &modelMIP, &ultimaColuna, geralSoluction, quantidade, &m);
		if (error) goto QUIT;

		//printLog(2);

	}

	//	printf("\nMontou modelo\n"); system("pause");

		//###############################################################################################################

		/*	error = GRBoptimize(model);
		if (error) { printf("\nErro ao otimizar o modelo RMP\n"); goto QUIT; }

		//printLog(3);

		//	printf("\nResolveu modelo\n");	system("pause");
		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, totalTipos, PiValue);
		if (error) { printf("\nErro ao obter os valores das variaveis duais\n"); goto QUIT; }

		//printLog(4);

		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
		if (error) { printf("\nErro ao inserir obter o valor da funcao objetivo - Primeira Iteracao\n"); goto QUIT; }

	//	printf("\nObteve os valores duais\n"); system("pause");

		idColuna = m - 1;

		if (classe == 9)qtdTiposBin = calloc(totalTipos, sizeof(int));

		nova_m = m;

		int variaveisArredondadas = 0;

		while (custo_reduzido_negativo == 1 && totalTime < timeLimitRMP) {

			//testeIntegridadeColunas(primeiraColuna, ultimaColuna, m);

			auxBin = CriaNovosPadroes(combinacoes[bestCombination][0], itens, heuristicas, ordem, PiValue, quantidade, tamanho, classe, nextStep);

			//printLog(5);

		//	printf("\nCriou novos padroes\n");system("pause");

			custo_reduzido_negativo = 0; //se foi encontrado custo reduzido negativo

			while (auxBin) {

				//printf("\nBin %d\n", auxBin->bin.idt);system("pause");

				reducedCost = 1.0; //valor do custo reduzido do padr�o atual

				auxPack = auxBin->bin.conteudo;

				//printLog(6);

				if (classe == 9) {

					zeraQtdTiposBin(qtdTiposBin, totalTipos);

					//printLog(7);

					totalTiposBin = 0;

				}

				while (auxPack && reducedCost >= 0.0) { //calcula custo reduzido

					if (classe == 9) {

						k = itens[auxPack->pack.id][5];

						//printLog(8);

						if (qtdTiposBin[k] == 0) totalTiposBin++;

						qtdTiposBin[k] = qtdTiposBin[k] + 1;

						//printLog(9);

					}
					else {

						k = auxPack->pack.id;

						//printLog(10);

					}

					if (PiValue[k] > 0.0) reducedCost = reducedCost - PiValue[k];

					//printLog(11);

					auxPack = auxPack->proximo;

					//printLog(12);
				}

				while (auxPack && classe == 9) {

					k = itens[auxPack->pack.id][5];

					//printLog(13);

					if (qtdTiposBin[k] == 0) totalTiposBin++;

					//printLog(14);

					qtdTiposBin[k] = qtdTiposBin[k] + 1;

					//printLog(15);

					auxPack = auxPack->proximo;

					//printLog(16);

				}

				if (reducedCost < 0.0) { //verifica se o custo reduzido � negativo

					//printLog(17);
					//printf("\nCusto-Reduzido Negativo\n", auxBin->bin.idt);	system("pause");

					//testeQtdTiposBin(qtdTiposBin, totalTipos, totalTiposBin, auxBin->bin.conteudo, itens);

					ultimaColuna = insereColuna(auxBin, ultimaColuna, idColuna);

					//printLog(18);

					//printf("\nInseriu coluna\n", auxBin->bin.idt);system("pause");

					if (classe == 9) {

						k = totalTiposBin;

						//printLog(19);

					}
					else {

						k = auxBin->bin.qtdeItens;

						//printLog(20);
					}

					alocaMemoriaUmaVariavel(&cind, &cval, k);

					if (classe == 9) {

						//printLog(21);

						nz = preencheUmaVariavelClass9(cind, cval, totalTipos, qtdTiposBin);

						//printLog(22);

						//if (nz != totalTiposBin) { printf("Erro! nz retorno funcao = %d, totaltiposBin %d", nz, totalTiposBin); system("pause"); }
						//else{printf("\nPassou no teste retorno nz preencheUmaVariavel\n");system("pause");}

					}
					else {

						//printLog(23);

						nz = preencheUmaVariavel(cind, cval, auxBin);

						//printLog(24);
					}

					error = GRBaddvar(model, nz, cind, cval, 1.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS, NULL);
					if (error) { printf("\nErro ao inserir nova coluna no RMP\n"); goto QUIT; }

					//printLog(25);
					//printf("\nAdicionou coluna ao RMP\n");system("pause");


					if(classe == 9){ error = GRBaddvar(modelMIP, nz, cind, cval, 1.0, 0.0, GRB_INFINITY, GRB_INTEGER, NULL); }
					else{ error = GRBaddvar(modelMIP, nz, cind, cval, 1.0, 0.0, GRB_INFINITY, GRB_BINARY, NULL); }

					if (error) { printf("\nErro ao inserir nova coluna no MIP\n"); goto QUIT; }

					//printLog(26);
					//printf("\nAdicionou coluna ao MIP\n");system("pause");

					nova_m = nova_m + 1; // aumenta uma vari�vel

					idColuna++;

					custo_reduzido_negativo = 1;

					auxBin = auxBin->proximo;

					//printLog(27);

				}
				else {

					//printLog(28);

					topoBin = auxBin;
					auxBin = auxBin->proximo;
					topoBin->proximo = NULL;
					freeMemoryBin(topoBin);

					//printLog(29);

				}

				//if (auxBin) { printf("\Vai passar para a proxima Bin\n"); system("pause"); }
				//else { printf("\Nao ha mais bin\n"); system("pause"); }

			}

			if (custo_reduzido_negativo == 0 || actualIteration > NIterations) {

				//printLog(30);

				custo_reduzido_negativo = 0;

		//		printf("\nVai arredondar uma variavel\n"); system("pause");

				if (variaveisArredondadas < (int)ceil(objval)) {

					//printLog(31);

					error = obtemResultado(model, &valorVariaveis, m);
					if (error) { printf("\nErro ao obter o resultado\n"); goto QUIT; }

					//printLog(32);

					//printf("\nObteve o resultado\n"); system("pause");

					error = fixaUmaVariavel(model, valorVariaveis, m, &custo_reduzido_negativo);
					if (error) goto QUIT;

					//printLog(33);

					if (custo_reduzido_negativo) variaveisArredondadas++;

					free(valorVariaveis);

				}

			}

			if (custo_reduzido_negativo == 1) {//encontrou-se padr�o com custo reduzido negativo

		//		printf("\nVai otimizar o modelo linear\n"); system("pause");

				//printLog(34);

				error = GRBoptimize(model);
				if (error) { printf("\nErro ao otimizar o RMP\n"); goto QUIT; }

				//printLog(35);

				error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);
				if (error) goto QUIT;

				//printLog(36);

				if (status != GRB_OPTIMAL) {

					fprintf(stderr, "Error: it isn't optimal\n");
					goto QUIT;
				}

		//		printf("\Otimizou o modelo\n"); system("pause");

				//printLog(37);

				GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, totalTipos, PiValue);
				if (error) { printf("\nErro ao obter os valores das variaveis duais (2)\n"); goto QUIT; }

				//printLog(38);

		//		printf("\Obteve o valor dos duais\n"); system("pause");

				//novo c�digo

				error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
				if (error) { printf("\nErro ao inserir obter o valor da funcao objetivo\n"); goto QUIT; }

				//printLog(39);

		//		printf("\Obteve o resultado\n"); system("pause");


				m = nova_m;

				if (referenceImprovement - objval >= ZMin) {

					//printLog(40);

					referenceImprovement = objval;

					actualIteration = 0;

					}
				else {

				//printLog(41);


				actualIteration++;

		//		if (actualIteration > NIterations) { printf("\nLimite de Iteracao sem melhoria\n"); system("pause"); }

						}
		//		printf("\Atualiza iteracao\n"); system("pause");

			}

			totalTime = time(NULL) - initialTime;

		//	if (totalTime > timeLimitRMP) { printf("\nTempo de execucao esgotado\n"); system("pause"); }

			//printLog(42);


		}

		m = nova_m;

		custo_reduzido_negativo = 1;

		tempoUltimaMelhoria = time(NULL) - initialTime;

		teta = 0.75;

		while (custo_reduzido_negativo == 1 && variaveisArredondadas < (int)ceil(objval)) {

			//printLog(43);

			custo_reduzido_negativo = 0;

			error = GRBoptimize(model);
			if (error) { printf("\nErro ao otimizar o RMP\n"); goto QUIT; }

			//printLog(44);

			error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);
			if (error) goto QUIT;

			//printLog(45);


			if (status != GRB_OPTIMAL) {

			fprintf(stderr, "Error: it isn't optimal\n");
			goto QUIT;
		}

		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
		if (error) { printf("\nErro ao obter o valor da funcao objetivo\n"); goto QUIT; }

		error = obtemResultado(model, &valorVariaveis, m);
		if (error) { printf("\nErro ao obter o resultado\n"); goto QUIT; }

		//printLog(46);

		error = fixaUmaVariavel(model, valorVariaveis, m, &custo_reduzido_negativo);
		if (error) goto QUIT;

		//printLog(47);

		if (custo_reduzido_negativo) variaveisArredondadas++;

		free(valorVariaveis);

	}


	//printLog(48);


	totalTime = time(NULL) - initialTime;

	roudingSoluction = empilharSoluction(roudingSoluction, NULL, 0, (int)ceil(objval), 0.0, totalTime);

	//printLog(49);

	//################################################################################################################

	*/

	totalTime = time(NULL) - initialTime;

	parcialTime = timeLimitMIP - totalTime;

	error = GRBsetdblparam(GRBgetenv(modelMIP), GRB_DBL_PAR_TIMELIMIT, (double)parcialTime);
	if (error) { printf("\nErro ao configurar o tempo limite MIP\n"); goto QUIT; }

	//printLog(50);


	//	printf("\nVai resolver o problema inteiro\n"); system("pause");

	error = GRBoptimize(modelMIP);
	if (error) { printf("\nErro ao otimizar o MIP\n"); goto QUIT; }

	//printLog(51);


	error = GRBgetintattr(modelMIP, GRB_INT_ATTR_STATUS, &status);

	if (status != GRB_OPTIMAL && status != GRB_TIME_LIMIT) {

		fprintf(stderr, "Error: it isn't optimal\n");
		goto QUIT;

	}

	error = GRBgetdblattr(modelMIP, GRB_DBL_ATTR_OBJVAL, &objval);
	if (error) { printf("\nErro ao obter o valor da funcao objetivo MIP\n"); goto QUIT; }

	error = obtemResultado(modelMIP, &valorVariaveis, m);
	if (error) { printf("\nErro ao obter resultado MIP\n"); goto QUIT; }

	//	}

//	printf("\nVai validar solucao\n"); system("pause");

	totalTime = time(NULL) - initialTime;

	//	totalTime = totalTime - roudingSoluction->soluction.time + tempoUltimaMelhoria;

	erro = validSoluction(bestOfBest->soluction.bins, tamanho, quantidade, bestOfBest->soluction.value);

	if (erro > 0) {

		printf("\n\tSolucao ACO nao gravada! Encontrado %i erros!\n", erro);

	}
	else {

		saveFile(bestOfBest->soluction.bins, classe, quantidade, instancia, bestOfBest->soluction.value, bestOfBest->soluction.time, semente);

	}

	salvaSolucao(&columnSoluction, primeiraColuna, valorVariaveis, m, (int)objval, totalTime);

	free(valorVariaveis);

	if (classe != 9) {


		excluiItensDuplicados(columnSoluction, itens, quantidade);

		erro = validSoluction(columnSoluction->soluction.bins, tamanho, quantidade, columnSoluction->soluction.value);


	}
	else {

		erro = validSoluctionClass9(columnSoluction->soluction.bins, tamanho, columnSoluction->soluction.value, totalTipos, qtdTipos, itens);

	}

	if (erro > 0) {

		printf("\n\tSolucao CG com ressalvas! Encontrado %i erros!\n", erro);

	}

	//saveFile(NULL, classe, quantidade, instancia, roudingSoluction->soluction.value, roudingSoluction->soluction.time, semente);

	saveFile(columnSoluction->soluction.bins, classe, quantidade, instancia, columnSoluction->soluction.value, columnSoluction->soluction.time, semente);

	printf("\n\tSolucao criada com sucesso! \n\tResultado registrado no arquivo BPP-0.1.6.csv\n");

	printf("\nFinalizado %d_%d_%d\n", classe, quantidade, instancia);

	printf("\n\n\t++++++++++++++++++++++++++++++++++++++++++++\n");


	//#####################################################################################


/*	if (classe == 9) {

		for (k = 0; k < totalTipos; k++) {

			free(qtdTipos[k]);

		}

		free(qtdTipos);

	}


	freeMemorySoluctionAll(columnSoluction);

	freeMemorySoluctionACO(geralSoluction);

	freeProbabilities(topoProbability);

	free(PiValue);
*/
	GRBfreemodel(modelMIP);

	/*	for (k = 0; k < quantidade; k++) {

			free(itens[k]);

		}

		for (k = 0; k < quantidadeOrdens; k++) {

			free(heuristicas[k]);

		}

		for (k = 0; k < quantidadeCombinacoes; k++) {

			free(pheromones[k]);

		}


		free(itens);
		free(ordem);
		free(bestOrder);
		free(nextStep);

		*/
		//Free environment
	GRBfreeenv(env);

QUIT:

	/* Error reporting */

	if (error) {

		printf("\nErro %d\n", error);

		printf("ERROR: %s\n", GRBgeterrormsg(env));
		system("pause");
		exit(1);
	}

	system("pause");

	return 0;

}
