#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "gurobi_c.h"

int  timeLimitACO, timeLimitCG, totalTimeLimit, NumberStart, method;
float QAnt, evaporation, QIncrease, pheromone, alfa, beta, nZao, QNZao, volume_Total;
int sr, spaceRule, orientationRule;
float teta;

typedef struct {

	double* mipStart;
	float value;

}BestStart;

typedef struct noBestStart {

	BestStart bestStart;
	struct noBestStart* proximo;

}NoBestStart;


BestStart preencheBestStart(int quantidade, float value) {

	BestStart bestStart;

	bestStart.mipStart = (double*)calloc(quantidade, sizeof(double));
	bestStart.value = value;

	return bestStart;

}

NoBestStart* insertBestStart(NoBestStart** atual, NoBestStart* topo, int quantidade, float value) {

	NoBestStart* novo = malloc(sizeof(NoBestStart));
	novo->bestStart = preencheBestStart(quantidade, value);

	*atual = novo;

	if (topo == NULL) {

		novo->proximo = NULL;
		return novo;

	}
	else if (topo->bestStart.value <= novo->bestStart.value) {

		novo->proximo = topo;
		return novo;

	}
	else {

		NoBestStart* anterior = topo, * aux = topo->proximo;

		while (aux && aux->bestStart.value > novo->bestStart.value) {

			anterior = aux;
			aux = aux->proximo;

		}

		anterior->proximo = novo;
		novo->proximo = aux;

		return topo;

	}

}

NoBestStart* deleteBestStart(NoBestStart* topo) {

	NoBestStart* aux = topo;
	topo = aux->proximo;
	free(aux->bestStart.mipStart);
	free(aux);
	return topo;
}


typedef struct {

	int id;
	float value;

}Probability;

typedef struct noProbability {

	Probability probability;
	struct noProbability* proximo;
	struct noProbability* permanentProximo;

}NoProbability;


Probability preencherProbability(int id) {

	Probability probability;

	probability.id = id;

	probability.value = 0.0;

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

void printProbability(NoProbability* topo) {

	NoProbability* aux = topo;

	printf("\n[");

	while (aux) {

		printf("%d: %f; ", aux->probability.id, aux->probability.value);

		aux = aux->permanentProximo;

	}

	printf("]\n");

}

void printProbabilityModified(NoProbability* topo) {

	NoProbability* aux = topo;

	printf("\n[");

	while (aux) {

		printf("%d: %f; ", aux->probability.id, aux->probability.value);

		aux = aux->proximo;

	}

	printf("]\n");

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
	int best;
	struct noBin* bins;

}Soluction;


typedef struct noSoluction {

	Soluction soluction;
	struct noSoluction* proximo; //ponteiro apontando para um estrutura do tipo noSoluction
	struct noSoluction* nextBestSoluction;

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

	printf("\nid: %i,  x: %i, y: %i, z: %i, orientation: %i, %i parallel to axis x, %i parallel to axis y, %i parallel to axis z,\n", pack.id, pack.x, pack.y, pack.z, pack.orient, pack.kx, pack.ky, pack.kz);
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
		inserir->space.y == aux->space.y && inserir->space.x < aux->space.x ||
		inserir->space.y == aux->space.y && inserir->space.x == aux->space.x && inserir->space.z < aux->space.z) {
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
		inserir->space.y == aux->space.y && inserir->space.z < aux->space.z ||
		inserir->space.y == aux->space.y && inserir->space.z == aux->space.z && inserir->space.x < aux->space.x) {
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
		inserir->space.x == aux->space.x && inserir->space.z < aux->space.z ||
		inserir->space.x == aux->space.x && inserir->space.z == aux->space.z && inserir->space.y < aux->space.y) {
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
		inserir->space.x == aux->space.x && inserir->space.y < aux->space.y ||
		inserir->space.x == aux->space.x && inserir->space.y == aux->space.y && inserir->space.z < aux->space.z) {
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
		inserir->space.z == aux->space.z && inserir->space.y < aux->space.y ||
		inserir->space.z == aux->space.z && inserir->space.y == aux->space.y && inserir->space.x < aux->space.x) {
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
		inserir->space.z == aux->space.z && inserir->space.x < aux->space.x ||
		inserir->space.z == aux->space.z && inserir->space.x == aux->space.x && inserir->space.y < aux->space.y) {
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

	}else {

		free(inserir);

	}

}

NoSpace* criarSpace(int x, int y, int z, int indexBin, int maxX, int maxY, int maxZ) {

	NoSpace* aux, * novo = malloc(sizeof(NoSpace));

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

//******************************procedimentos Soluction*****************************

Soluction preencherSoluction(int value, float utilization, int time, int id, NoBin* listaBin) {


	Soluction soluction;

	soluction.id = id;
	soluction.value = value;
	soluction.utilization = utilization;
	soluction.time = time;
	soluction.bins = listaBin;
	soluction.best = 0;


	return soluction;
}

NoSoluction* empilharSoluction(NoSoluction* listaSoluction, NoBin* listaBin, int id, int value, float utilization, int time) {
	NoSoluction* novo = malloc(sizeof(NoSoluction));

	if (novo) {

		novo->soluction = preencherSoluction(value, utilization, time, id, listaBin);//salva o dado pessoa na struct novo
		novo->proximo = listaSoluction; //salva ponteiro próximo como o último topo
		novo->nextBestSoluction = NULL;


		return novo;

	}
	else {

		printf("Não foi possível alocar memória");
		return NULL;
	}

}

NoSoluction* empilharBestSoluction(NoSoluction* listaBestSoluction, NoSoluction* novo) {

	novo->nextBestSoluction = listaBestSoluction; //salva ponteiro próximo como o último topo

	return novo;

}


NoSoluction* findSoluction(int id, NoSoluction* lista) { //encontra o ponteiro para um compartimento (Bin) a partir de seu index

	NoSoluction* usedSoluction = lista;

	while (usedSoluction && (*usedSoluction).soluction.id != id) {

		usedSoluction = usedSoluction->proximo;

	}


	return usedSoluction;

}


void freeMemorySoluction(NoSoluction* geralSoluction) {

	freeMemoryBin(geralSoluction->soluction.bins);

	free(geralSoluction);

}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//este grupo de  funções ordena os itens de acordo com o cálculo de próximo passo e retorna um vetor com os índices ordenados

void ordenaValores(float** nextStep, int* ordem, int quantidade) {


	void quicksort(int inicio, int fim, float** nextStep, int* ordem);
	int particao(int inicio, int fim, float** nextStep, int* ordem);
	int contador;

	for (contador = 0; contador < quantidade; contador++) {

		ordem[contador] = contador;

	}


	quicksort(0, quantidade - 1, nextStep, ordem);

}


void quicksort(int inicio, int fim, float** nextStep, int* ordem) {

	if (inicio < fim) {

		int pivo = particao(inicio, fim, nextStep, ordem);
		quicksort(inicio, pivo - 1, nextStep, ordem);
		quicksort(pivo + 1, fim, nextStep, ordem);

	}


}

int particao(int inicio, int fim, float** nextStep, int* ordem) {
	int i = inicio;
	int j;
	int temporario;

	for (j = inicio; j < fim; j++) {

		if (nextStep[ordem[j]][0] >= nextStep[ordem[fim]][0]) {

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


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//importa a instância do arquivo txt e salva no tabela itens
void importaArquivo(int** itens, int classe, int quantidade, int instancia, int tamanho[]) {

	volume_Total = 0.0;

	char texto1[15];
	char texto2[] = "_";
	char texto3[5];
	char texto4[3];
	char texto5[] = ".txt";

	sprintf(texto1, "%d", classe);
	sprintf(texto3, "%d", quantidade);
	sprintf(texto4, "%d", instancia);


	strcat(texto1, texto2);
	strcat(texto1, texto3);
	strcat(texto1, texto2);
	strcat(texto1, texto4);
	strcat(texto1, texto5);

	printf("\n\n\t++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("\n\tInstancia utilizada %s\n\n", texto1);

	int contador;
	int maiorVolume = 0;

	FILE* file;
	file = fopen(texto1, "r");

	if (classe == 9) {

		fscanf(file, "%i", &tamanho[0]);
		fscanf(file, "%i", &tamanho[1]);
		fscanf(file, "%i", &tamanho[2]);

	}

	for (contador = 0; contador < quantidade; contador++) {

		fscanf(file, "%i", &itens[contador][0]);
		fscanf(file, "%i", &itens[contador][1]);
		fscanf(file, "%i", &itens[contador][2]);
		fscanf(file, "%i", &itens[contador][3]);

		if (classe > 8) {


			if (maiorVolume < itens[contador][1] * itens[contador][2] * itens[contador][3]) {

				maiorVolume = itens[contador][1] * itens[contador][2] * itens[contador][3];

			}

			if (classe == 10) volume_Total = volume_Total + (((float)itens[contador][1]) / 100.0 * ((float)itens[contador][2]) / 100.0 * ((float)itens[contador][3]) / 100.0);

			if (classe == 9) fscanf(file, "%i", &itens[contador][5]);

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

definePointAxis(int pointAxis[], int packed, int x, int y, int z) {

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


	NoSpace* atual=space, *anterior = space->anteriorGeral;

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

	NoSpace* deletar, * anterior, * aux = usedBin->bin.spaces;

	while (aux) {

		resolveOverlap(aux, item);

		if (aux->space.volume <= 0) {

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

	if (aux && aux->space.volume <= 0) {

		free(aux);
		*spaceXLeft = NULL;

	}

	aux = *spaceXDown;

	if (aux && aux->space.volume <= 0) {

		free(aux);
		*spaceXDown = NULL;

	}

	aux = *spaceYDown;

	if (aux && aux->space.volume <= 0) {

		free(aux);
		*spaceYDown = NULL;

	}

	aux = *spaceYBack;

	if (aux && aux->space.volume <= 0) {

		free(aux);
		*spaceYBack = NULL;

	}

	aux = *spaceZLeft;

	if (aux && aux->space.volume <= 0) {

		free(aux);
		*spaceZLeft = NULL;

	}

	aux = *spaceZBack;

	if (aux && aux->space.volume <= 0) {

		free(aux);
		*spaceZBack = NULL;

	}

}


NoBin* binPack(int** itens, int* ordem, int quantidade, int tamanho[]) {

	int i, o, x, y, z, id, packed = 0, volumes;
	int idBin;
	int pointAxis[3];
	NoBin* topoBin, * usedBin, * bestBins = NULL;
	NoSpace* aux, * topoSpaceGeral;
	NoPack* topoPack;
	NoSpace* spaceXLeft, * spaceXDown, * spaceYDown, * spaceYBack, * spaceZLeft, * spaceZBack, * novoSpace;

	for (o = -1; o < orientationRule; o++) {

		for (sr = 0; sr < spaceRule; sr++) {

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

				}

			}

			freeMemorySpace(topoSpaceGeral);

			if (bestBins == NULL || topoBin->bin.idt < bestBins->bin.idt) {

				freeMemoryBin(bestBins);
				bestBins = topoBin;

			}
			else {

				freeMemoryBin(topoBin);

			}

		}


	}

	return bestBins;

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



//valida solução final para verificar todos os itens estão no limite da Bin e se não há overleap entre nenhum item
int validSoluction(NoBin* listaBin, int tam[], int quantidade, int solucao) {

	int erro = 0, totalItens = 0, totalBin = 0, i;
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

				printf("\n\tItem %d duplicated!\n");

				erro++;


			}


			totalItens++;

			if ((*auxPack).pack.x + (*auxPack).pack.kx > tam[0] || (*auxPack).pack.y + (*auxPack).pack.ky > tam[1] || (*auxPack).pack.z + (*auxPack).pack.kz > tam[2]) {

				printf("\n\tInvalid Solution! Item item exceeds the bin limit!\n");
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

					printf("\n\tInvalid Solution! Itens overleap!\n");
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

		printf("\n\tItems loaded %d! Items quantity %d!\n", totalItens, quantidade);

		erro++;

	}

	if (totalBin != solucao) {

		printf("\n\tBins used %d! Solution %d!\n", totalBin, solucao);

		erro++;

	}

	free(used);

	return erro;

}


//salva o arquivo com os resultados
void saveFile(NoBin* listaBin, int classe, int quantidade, int instancia, int objetivo, int tempo, unsigned int semente) {

	FILE* pont_arqu;
	NoBin* auxBin = listaBin;
	NoPack* auxPack;
	NoPack* auxPack1;
	float wasted = 0;


	pont_arqu = fopen("results.csv", "a");


	if (classe == 10) {

		wasted = ((float)objetivo) * 36.24556364 - volume_Total;

	}

	fprintf(pont_arqu, "%d;%d;%d;%d;%d;%d;%.2f; [", semente, classe, quantidade, instancia, objetivo, tempo, wasted);

	while (auxBin) {// percorre todas as bins utilizadas

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

	fclose(pont_arqu);

}



//**********************************Conjunto de funcoes do Algoritmo Colônia de Formigas*************************************

void printParams() {

	printf("\n\tValores dos parametros");
	printf("\n\tQAnt\t\t\t %f", QAnt);
	printf("\n\talpha value\t\t\t %f", alfa);
	printf("\n\tbeta value\t\t\t %f", beta);
	printf("\n\tevaporation proportion\t\t %f", evaporation);
	printf("\n\taddition proportion\t\t %f", QIncrease);
	printf("\n\tpheromone value\t\t %f", pheromone);
	printf("\n\tQNZao\t\t\t %f", QNZao);
	printf("\n\tACO time limit\t %d", timeLimitACO);
	printf("\n\t RMP time limit\t %d", timeLimitCG);
	printf("\n\tTotal time limit\t %d", totalTimeLimit);
	printf("\n\tMIP Starts quantity\t %d", NumberStart);
	printf("\n\tResidual Space rules quantity\t%d", spaceRule);
	printf("\n\tOrientation rules quantity\t%d", orientationRule);
	printf("\n\tTheta value (fixed CG variables)\t%f", teta);


}

void obtemValoresParametros() {

	int opcao = 11;

	while (opcao > 0) {

		printf("\nType:\n");
		printf("\t0 - To return.\n");
		printf("\t1 - To change QAnt.\n");
		printf("\t2 - To change Alpha.\n");
		printf("\t3 - To change Beta.\n");
		printf("\t4 - To change evaporation proportion.\n");
		printf("\t5 - To change addition proportion.\n");
		printf("\t6 - To change pheromone value.\n");
		printf("\t7 - To change ACO time limit.\n");
		printf("\t8 - To change RMP time limit (CG).\n");
		printf("\t9 - To change total time limit.\n");
		printf("\t10 - To change MIP Start quantity.\n");
		printf("\t11 - To change residual space rule quantity.\n");
		printf("\t12 - To change orientation rule quantity.\n");
		printf("\t13 - To change theta value.\n");


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
			printf("\n\tType ACO time limit in seconds(Ex: 10)\n\t");
			scanf("%d", &timeLimitACO);
			break;

		case 8:
			printf("\n\tType RMP (CG) time limit in seconds (Ex: 10)\n\t");
			scanf("%d", &timeLimitCG);
			break;

		case 9:
			printf("\n\tType total time limit in seconds (Ex: 10)\n\t");
			scanf("%d", &totalTimeLimit);
			break;

		case 10:
			printf("\n\tType MIP Stars quantity (Ex: 5)\n\t");
			scanf("%d", &NumberStart);
			break;

		case 11:
			printf("\n\tType residual space rules quantity (Ex: 1-smaller volume, from 2 to 7 - smaller volume and extreme-point location)\n\t");
			scanf("%d", &spaceRule);
			break;

		case 12:
			printf("\n\tType orientation rules quantity (Ex: Ex: 0-Smaller difference, 1 a 6 - smaller difference and preference sequence)\n\t");
			scanf("%d", &orientationRule);
			break;

		case 13:
			printf("\n\tType theta value (Ex: 0.7)\n\t");
			scanf("%f", &teta);
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
	file = fopen("parametros.txt", "r");

	fscanf(file, "%f", &QAnt);
	fscanf(file, "%f", &alfa);
	fscanf(file, "%f", &beta);
	fscanf(file, "%f", &evaporation);
	fscanf(file, "%f", &QIncrease);
	fscanf(file, "%f", &pheromone);
	fscanf(file, "%f", &QNZao);
	fscanf(file, "%d", &timeLimitACO);
	fscanf(file, "%d", &timeLimitCG);
	fscanf(file, "%d", &totalTimeLimit);
	fscanf(file, "%d", &NumberStart);
	fscanf(file, "%d", &spaceRule);
	fscanf(file, "%d", &orientationRule);
	fscanf(file, "%f", &teta);

	fclose(file);

}


void updateParams() {

	FILE* pont_arqu;

	pont_arqu = fopen("parametros.txt", "w");

	//	fprintf(pont_arqu, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%lf\n", QAnt, alfa, beta, evaporation, QIncrease, pheromone, QNZao, timeLimit, mipTimeLimit);

	fprintf(pont_arqu, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", QAnt, alfa, beta, evaporation, QIncrease, pheromone, QNZao, timeLimitACO, timeLimitCG, totalTimeLimit, NumberStart, spaceRule, orientationRule, teta);

	fclose(pont_arqu);

	system("cls");

	printParams();

}

//calcula a heurística. Precisa ser calculado uma única vez
float calcHeuristic(int** itens, float** nextStep, int quantidade) {

	int j;
	float totalVolume = 0.0;

	for (j = 0; j < quantidade; j++) {


		totalVolume = totalVolume + ((float)(itens[j][1] * itens[j][2] * itens[j][3]));
		nextStep[j][2] = pow((float)(itens[j][1] * itens[j][2] * itens[j][3]) / nZao, beta);
		nextStep[j][3] = pheromone;

	}

	return totalVolume;

}


//calcula o denominador. Calculado a cada iteração
float calcDenominator(float** nextStep, int quantidade, NoProbability* topo) {

	NoProbability* aux = topo;
	int j;
	float denominator = 0;

	for (j = 0; j < quantidade; j++) {

		nextStep[j][1] = pow(nextStep[j][3], alfa);

		aux->probability.value = nextStep[j][1] * nextStep[j][2];

		denominator = denominator + aux->probability.value;

		aux = aux->permanentProximo;

	}


	return denominator;

}


//calcula a probabilidade, deve ser calculado para cada formiga
void calcProbality(int quantidade, float denominator, NoProbability* topo, int* ordem) {

	int j;
	float luck, acumulated;
	NoProbability* aux, * anterior, * novotopo = topo;
	float newdenominator = denominator;

	for (j = 0; j < quantidade; j++) {

		luck = ((float)rand() / (float)RAND_MAX);

		aux = novotopo;

		anterior = NULL;

		acumulated = aux->probability.value / newdenominator;

		while (aux->proximo && acumulated < luck) {

			anterior = aux;

			aux = aux->proximo;

			acumulated = acumulated + (aux->probability.value / newdenominator);

		}

		ordem[j] = aux->probability.id;

		newdenominator = newdenominator - aux->probability.value;

		deleteProbability(&novotopo, anterior, aux);

	}

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


void copyVector(int* ordem, int* bestOrder, int quantidade) {

	int i;

	for (i = 0; i < quantidade; i++) {

		bestOrder[i] = ordem[i];

	}
}

void updatePheromone(int* bestOrder, float** nextStep, int quantidade, float increase, float decrease) {

	int i;

	for (i = 0; i < quantidade; i++) {


		nextStep[bestOrder[i]][3] = (nextStep[bestOrder[i]][3] * (((float)1) - evaporation)) + pheromone - (((float)i) * decrease) + ((float)increase) * QIncrease;

	}


}

NoProbability* montaPilhaProbabilidade(int quantidade) {

	NoProbability* aux = NULL;
	int i;

	for (i = quantidade - 1; i >= 0; i--) {

		aux = insereProbability(aux, i);

	}

	return aux;

}


void alocaMemoriaRestricao(double** rhs, char** sense, char** senseMIP, int n) {

	*rhs = (double*)malloc(n * sizeof(double)); //reserva memória para o lado direito da restrição (double)

	*sense = (char*)malloc(n * sizeof(char)); //reserva memória para o sinal das restrições (char)

	*senseMIP = (char*)malloc(n * sizeof(char)); //reserva memória para o sinal das restrições (char)

}

void criaRestricoes(double* rhs, char* sense, char* senseMIP, int n) {

	int i;

	for (i = 0; i < n; i++) {

		rhs[i] = 1.0; //armazena o lado direito da restrição

		sense[i] = GRB_EQUAL; //armazena o sinal da restrição

		senseMIP[i] = GRB_GREATER_EQUAL;

	}

}

void alocaMemoriaVariaveis(int** cbeg, int** clen, int** cind, double** cval, double** lb, double** obj, char** ctype, char** ctypeMIP, int m, int nz) {

	*cbeg = (int*)malloc(m * sizeof(int)); //armazena a quantidade de variáveis (inteiro)

	*clen = (int*)malloc(m * sizeof(int)); //armazena a quantidade variáveis na restrição (inteiro)

	*cind = (int*)malloc(nz * sizeof(int)); //armazena os índices da restrição (inteiro)

	*cval = (double*)malloc(nz * sizeof(double)); //armazena os valores dos multiplicadores das variáveis (double)

	*lb = (double*)malloc(m * sizeof(double)); //armazena o limite inferior das variáveis (double)

	*obj = (double*)malloc(m * sizeof(double)); //armazena o multiplicador das variáveis na função objetivo (double)

	*ctype = (char*)malloc(m * sizeof(char)); //armazena o tipo de variaval (Integer, Binary, Continuous)

	*ctypeMIP = (char*)malloc(m * sizeof(char)); //armazena o tipo de variaval (Integer, Binary, Continuous)

}

void alocaMemoriaUmaVariavel(int** cind, double** cval, NoBin* auxBin) {

	*cind = (int*)malloc(auxBin->bin.qtdeItens * sizeof(int)); //armazena os índices da restrição (inteiro)

	*cval = (double*)malloc(auxBin->bin.qtdeItens * sizeof(double)); //armazena os valores dos multiplicadores das variáveis (double)

}


void criaFuncaoObjetivo(double* lb, double* obj, int m) {

	int j;

	for (j = 0; j < m; j++) {
		//for x variable

		lb[j] = 0.0; //lower bound da variável x

		obj[j] = 1.0; //valor do objetivo = c

	}

}

NoBestStart* preencheVariaveis(NoBin** ultimaColuna, NoSoluction* soluctions, int* cbeg, int* cind, double* cval, int* clen, char* ctype, char* ctypeMIP, int quantidade) {

	int nz = 0;

	int numeroStart = 0, i = 0;

	NoBestStart* bestAtual, * topoBestStart = NULL;

	NoSoluction* auxSoluction = soluctions;

	NoBin* auxBin, * binAnterior = NULL;

	NoPack* auxPack;

	while (auxSoluction) {

		if (numeroStart < NumberStart) {

			topoBestStart = insertBestStart(&bestAtual, topoBestStart, quantidade, auxSoluction->soluction.utilization);
			numeroStart++;


		}
		else if (auxSoluction->soluction.utilization < topoBestStart->bestStart.value) {

			topoBestStart = insertBestStart(&bestAtual, topoBestStart, quantidade, auxSoluction->soluction.utilization);
			topoBestStart = deleteBestStart(topoBestStart);


		}
		else {

			bestAtual = NULL;
		}

		auxBin = auxSoluction->soluction.bins;

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

			if (bestAtual) {//marca a variavel no MIP Start

				bestAtual->bestStart.mipStart[i] = 1.0;

			}


			ctype[i] = GRB_CONTINUOUS;

			ctypeMIP[i] = GRB_BINARY;

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

	//testar esse procedimento
	*ultimaColuna = binAnterior;

	return topoBestStart;

}


int preencheUmaVariavel(int* cind, double* cval, NoBin* auxBin) {

	int nz = 0;

	NoPack* auxPack = auxBin->bin.conteudo;

	while (auxPack) {

		//testeMemoriaVariavel(nz, auxBin->bin.qtdeItens);

		cind[nz] = auxPack->pack.id;

		cval[nz] = 1.0;

		nz++;

		auxPack = auxPack->proximo;
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
	char* ctype = NULL; //tipo da variável
	char* ctypeMIP = NULL; //tipo da variável

	NoBestStart* topoBestStart, * auxBestStart;

	NoSoluction* auxSoluction = geralSoluction;
	NoBin* topoBin = NULL;

	alocaMemoriaRestricao(&rhs, &sense, &senseMIP, n);

	criaRestricoes(rhs, sense, senseMIP, n);

	while (auxSoluction) {

		*m = *m + auxSoluction->soluction.value; // quantidade de variáveis

		nz = nz + n; //quantidade de indices nao zero nas restricoes

		auxSoluction = auxSoluction->proximo;

	}

	alocaMemoriaVariaveis(&cbeg, &clen, &cind, &cval, &lb, &obj, &ctype, &ctypeMIP, *m, nz);

	criaFuncaoObjetivo(lb, obj, *m);

	topoBestStart = preencheVariaveis(ultimaColuna, geralSoluction, cbeg, cind, cval, clen, ctype, ctypeMIP, *m);

	error = GRBsetintparam(env, "StartNumber", -1);
	if (error) return error;

	error = GRBloadmodel(env, model, "Linear", *m, n,
		GRB_MINIMIZE, 0.0, obj, sense, rhs,
		cbeg, clen, cind, cval, lb, NULL,
		ctype, NULL, NULL);
	if (error) return error;

	error = GRBloadmodel(env, modelMIP, "MIP", *m, n,
		GRB_MINIMIZE, 0.0, obj, senseMIP, rhs,
		cbeg, clen, cind, cval, lb, NULL,
		ctypeMIP, NULL, NULL);
	if (error) return error;

	error = GRBupdatemodel(*modelMIP);
	if (error) return error;

	auxBestStart = topoBestStart;

	while (auxBestStart) {

		error = GRBsetdblattrarray(*modelMIP, "Start", 0, *m, auxBestStart->bestStart.mipStart);
		if (error) return error;

		topoBestStart = auxBestStart;
		auxBestStart = auxBestStart->proximo;
		free(topoBestStart);


	}

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

int obtemResultado(GRBmodel* model, double** valorVariaveis, int m) {

	int erro;

	*valorVariaveis = (double*)malloc(m * sizeof(double));

	erro = GRBgetdblattrarray(model, "X", 0, m, *valorVariaveis);

	return erro;
}


void salvaSolucao(NoSoluction** columnSoluction, NoBin* primeiraColuna, double* valorVariaveis, int m, int valor, int totalTime) {

	NoBin* auxBin = primeiraColuna, * topoBin = NULL;

	int j;

	for (j = 0; j < m; j++) {

		if (valorVariaveis[j] > 0.0) {

			auxBin->proximo = topoBin;
			topoBin = auxBin;

		}

		auxBin = auxBin->proximoColuna;

	}

	*columnSoluction = empilharSoluction(*columnSoluction, topoBin, 0, valor, 0.0, totalTime);

}


void excluiItensDuplicados(int* usedItens, NoSoluction* columnSoluction, int** itens) {

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


}

void  valoresIniciaisDeOrdem(double* PiValue, float** nextStep, int* ordem, int quantidade) {

	int j;

	for (j = 0; j < quantidade; j++) {

		ordem[j] = j;
		nextStep[j][0] = (float)PiValue[j];
		nextStep[j][1] = (float)PiValue[j];

	}

}

NoBin* CriaNovosPadroes(int** itens, float** nextStep, int* ordem, double* PiValue, int quantidade, int tamanho[], int totalTime) {

	NoBin* auxBin = NULL, * auxBin1 = NULL;
	NoBin* topoBin = NULL;
	int k;

	//Salva no vetor Ordem Apenas os itens com valores Duais Positivos
	valoresIniciaisDeOrdem(PiValue, nextStep, ordem, quantidade);

	//ordena pelos valores das variáveis duais
	quicksort(0, quantidade - 1, nextStep, ordem);

	topoBin = binPack(itens, ordem, quantidade, tamanho);

	for (k = 0; k < quantidade; k++) {

		nextStep[ordem[k]][0] = nextStep[ordem[k]][1] / (float)(itens[ordem[k]][1] * itens[ordem[k]][2] * itens[ordem[k]][3]);

	}

	//ordena pelos valores das variáveis duais dividido pelo volume
	quicksort(0, quantidade - 1, nextStep, ordem);

	auxBin = binPack(itens, ordem, quantidade, tamanho);

	auxBin1 = auxBin;

	while (auxBin1->proximo) {

		auxBin1 = auxBin1->proximo;

	}

	auxBin1->proximo = topoBin;
	topoBin = auxBin;

	for (k = 0; k < quantidade; k++) {

		nextStep[ordem[k]][0] = nextStep[ordem[k]][1] / (float)(itens[ordem[k]][1] + itens[ordem[k]][2] + itens[ordem[k]][3]);

	}

	//ordena pelos valores das variáveis duais dividido pela soma dos lados
	quicksort(0, quantidade - 1, nextStep, ordem);

	auxBin = binPack(itens, ordem, quantidade, tamanho);

	auxBin1 = auxBin;

	while (auxBin1->proximo) {

		auxBin1 = auxBin1->proximo;

	}

	auxBin1->proximo = topoBin;
	topoBin = auxBin;

	for (k = 0; k < quantidade; k++) {

		nextStep[ordem[k]][0] = nextStep[ordem[k]][1] / (float)itens[ordem[k]][1];

	}

	//ordena pelos valores das variáveis duais dividido pelo lado x
	quicksort(0, quantidade - 1, nextStep, ordem);

	auxBin = binPack(itens, ordem, quantidade, tamanho);

	auxBin1 = auxBin;

	while (auxBin1->proximo) {

		auxBin1 = auxBin1->proximo;

	}

	auxBin1->proximo = topoBin;
	topoBin = auxBin;

	for (k = 0; k < quantidade; k++) {

		nextStep[ordem[k]][0] = nextStep[ordem[k]][1] / (float)itens[ordem[k]][2];

	}

	//ordena pelos valores das variáveis duais dividido pelo lado y
	quicksort(0, quantidade - 1, nextStep, ordem);

	auxBin = binPack(itens, ordem, quantidade, tamanho);

	auxBin1 = auxBin;

	while (auxBin1->proximo) {

		auxBin1 = auxBin1->proximo;

	}

	auxBin1->proximo = topoBin;
	topoBin = auxBin;

	for (k = 0; k < quantidade; k++) {

		nextStep[ordem[k]][0] = nextStep[ordem[k]][1] / (float)itens[ordem[k]][3];

	}

	//ordena pelos valores das variáveis duais dividido pelo lado z
	quicksort(0, quantidade - 1, nextStep, ordem);

	auxBin = binPack(itens, ordem, quantidade, tamanho);

	auxBin1 = auxBin;

	while (auxBin1->proximo) {

		auxBin1 = auxBin1->proximo;

	}

	auxBin1->proximo = topoBin;
	topoBin = auxBin;


	return topoBin;

}


int fixaUmaVariavel(GRBmodel* model, double* valorVariaveis, int m) {

	int k, id = -1, error = 0;
	double maiorValor = 0.0;

	for (k = 0; k < m; k++) {

		if (valorVariaveis[k] > (double)teta && valorVariaveis[k] > maiorValor) {

			id = k;
			maiorValor = valorVariaveis[k];

		}


	}

	if (id > -1) {

		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, id, 1.0);

	}

	free(valorVariaveis);

	return error;

}


int main() {

	int idSoluction, tempo, opcao, tamanho[3], classe, quantidade, instancia, totalTime, ants, parcialTime, actualAnt, erro;
	float lowerBound, denominator, decrease, utilization;
	NoBin* topoBin, * auxBin, * ultimaColuna = NULL, * primeiraColuna = NULL;
	NoSoluction* bestIteration, * bestOfBest, * geralSoluction;
	NoSoluction* columnSoluction = NULL;
	NoPack* auxPack;
	//Variaveis CG

	GRBenv* env = NULL; //criação do ambiente
	GRBmodel* model = NULL;
	GRBmodel* modelMIP = NULL; //criação do modelo
	int idColuna = 0; //id da última coluna adicionada
	int* cind = NULL;//contém os indices da variável na restrição
	double* cval = NULL;//indica o multiplicador da variavel na restrição
	double* valorVariaveis = NULL;

	int  nz = 0, i, j; //nz não zeros, i = x, j = contador
	int status;
	double objval;

	int error = 0, m = 0, n;

	int custo_reduzido_negativo = 1; // se há custo reduzido negativo
	int integer = 0; //se a solução é inteira
	double reducedCost; //valor do custo reduzido do padrão atual

	unsigned int semente = (unsigned)time(NULL);

	int** itens;
	int* ordem, * bestOrder, * usedItens;
	float** nextStep;
	double* PiValue;


	importaParametros();

	printf("\tVersion 0.0.12\t\tMethods ACO+CG\n\tDetails:");
	printf("\n\tMultiple Orientations -included Mahvash criteria- and multiple Residual Spaces\n");
	printf("\tIntegral Columns, Salve seed\n");
	printf("\tLarge Instances\n");
	printf("\n\tAuthor: Daniel Bento Maia\tDate: 03/08/2024\n");

	printParams();

	printf("\n\n\tTo change the parameters type 2; to solve one problem, type 1; to exit, type 0\n\t");
	scanf("%d", &opcao);

	if (opcao == 2) {

		obtemValoresParametros();

		updateParams();

		opcao = 1;
	}


	while (opcao == 1) {

		printf("\n\n\tPlease type class, quantity of itens and instance:");
		printf("\n\tFor example, to solve file 1_50_3.txt type 1 50 3\n\n\t");
		scanf("%d %d %d", &classe, &quantidade, &instancia);

		itens = malloc(quantidade * sizeof(int*));
		nextStep = malloc(quantidade * sizeof(float*));
		ordem = (int*)malloc(quantidade * sizeof(int));
		bestOrder = (int*)malloc(quantidade * sizeof(int));
		PiValue = (double*)malloc(quantidade * sizeof(double));
		usedItens = (int*)calloc(quantidade, sizeof(int));

		int k;

		for (k = 0; k < quantidade; k++) {

			itens[k] = malloc(6 * sizeof(int*));
			nextStep[k] = malloc(4 * sizeof(float*));

		}

		NoProbability* topoProbability = NULL;

		topoProbability = montaPilhaProbabilidade(quantidade);

		tamanhoInstancia(classe, tamanho);

		ants = ((float)quantidade) * QAnt;

		decrease = (float)pheromone / (float)quantidade;

		idSoluction = 0;
		totalTime = 0;
		topoBin = NULL;
		bestOfBest = NULL;
		geralSoluction = NULL;

		importaArquivo(itens, classe, quantidade, instancia, tamanho);

		error = GRBloadenv(&env, NULL);
		if (error) goto QUIT;

		srand(semente);

		tempo = time(NULL);

		firstOrientation(itens, quantidade);

		lowerBound = calcHeuristic(itens, nextStep, quantidade);

		lowerBound = lowerBound / ((float)(tamanho[0] * tamanho[1] * tamanho[2]));

		while (totalTime < timeLimitACO) {

			actualAnt = 0;

			denominator = calcDenominator(nextStep, quantidade, topoProbability);

			bestIteration = NULL;

			while (actualAnt <= ants) {

				calcProbality(quantidade, denominator, topoProbability, ordem);

				topoBin = binPack(itens, ordem, quantidade, tamanho);

				totalTime = time(NULL) - tempo;

				utilization = calcSmallerUtilization(topoBin, tamanho);

				if (bestIteration == NULL || bestIteration->soluction.utilization > (float)topoBin->bin.idt + (float)1 + utilization) {

					if (bestIteration != NULL) {


						freeMemorySoluction(bestIteration);
						bestIteration = NULL;


					}

					bestIteration = empilharSoluction(bestIteration, topoBin, idSoluction, topoBin->bin.idt + 1, (float)topoBin->bin.idt + (float)1 + utilization, totalTime);

					copyVector(ordem, bestOrder, quantidade);

				}
				else {

					freeMemoryBin(topoBin);


				}

				actualAnt++;

			}

			bestIteration->proximo = geralSoluction;

			geralSoluction = bestIteration;

			updatePheromone(bestOrder, nextStep, quantidade, lowerBound / bestIteration->soluction.utilization, decrease);

			if (bestOfBest == NULL || bestIteration->soluction.value < bestOfBest->soluction.value) {

				if (bestOfBest) {

					bestOfBest->soluction.best = 0;

				}

				bestOfBest = empilharBestSoluction(bestOfBest, bestIteration);

				bestOfBest->soluction.best = 1;

			}

		}

		primeiraColuna = geralSoluction->soluction.bins;

		error = montaModelo(env, &model, &modelMIP, &ultimaColuna, geralSoluction, quantidade, &m);
		if (error) goto QUIT;

		error = GRBoptimize(model);
		if (error) goto QUIT;

		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, quantidade, PiValue);
		if (error) goto QUIT;

		idColuna = m - 1;

		while (totalTime < timeLimitCG && custo_reduzido_negativo == 1) {

			auxBin = CriaNovosPadroes(itens, nextStep, ordem, PiValue, quantidade, tamanho, totalTime);

			custo_reduzido_negativo = 0; //se foi encontrado custo reduzido negativo

			while (auxBin) {

				reducedCost = 1.0; //valor do custo reduzido do padrão atual

				auxPack = auxBin->bin.conteudo;

				while (auxPack && reducedCost >= 0.0) { //calcula custo reduzido

					if (PiValue[auxPack->pack.id] > 0.0) reducedCost = reducedCost - PiValue[auxPack->pack.id];

					auxPack = auxPack->proximo;

				}

				if (reducedCost < 0.0) { //verifica se o custo reduzido é negativo

					ultimaColuna = insereColuna(auxBin, ultimaColuna, idColuna);

					alocaMemoriaUmaVariavel(&cind, &cval, auxBin);

					nz = preencheUmaVariavel(cind, cval, auxBin);

					error = GRBaddvar(model, nz, cind, cval, 1.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS, NULL);
					if (error) goto QUIT;

					//#######################################################################################################
					error = GRBaddvar(modelMIP, nz, cind, cval, 1.0, 0.0, GRB_INFINITY, GRB_BINARY, NULL);
					if (error) goto QUIT;

					m = m + 1; // aumenta uma variável

					idColuna++;

					custo_reduzido_negativo = 1;

					auxBin = auxBin->proximo;

				}
				else {

					topoBin = auxBin;
					auxBin = auxBin->proximo;
					topoBin->proximo = NULL;
					freeMemoryBin(topoBin);

				}


			}

			if (custo_reduzido_negativo == 1) {//encontrou-se padrão com custo reduzido negativo


				error = GRBoptimize(model);
				if (error) goto QUIT;

				GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, quantidade, PiValue);
				if (error) goto QUIT;

				//novo código
				error = obtemResultado(model, &valorVariaveis, m);
				if (error) goto QUIT;

				error = fixaUmaVariavel(model, valorVariaveis, m);
				if (error) goto QUIT;
				//novo código


			}

			totalTime = time(NULL) - tempo;

		}

		if (verificaSeEhInteiro(model, m)) {

			error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
			if (error) goto QUIT;

			error = obtemResultado(model, &valorVariaveis, m);
			if (error) goto QUIT;

		}
		else {

			//#################################################

			parcialTime = totalTimeLimit - (time(NULL) - tempo);

			error = GRBsetdblparam(GRBgetenv(modelMIP), GRB_DBL_PAR_TIMELIMIT, (double)parcialTime);
			if (error) goto QUIT;

			//Optimize model

			error = GRBoptimize(modelMIP);
			if (error) goto QUIT;

			totalTime = time(NULL) - tempo;

			error = GRBgetintattr(modelMIP, GRB_INT_ATTR_STATUS, &status);
			if (error) goto QUIT;

			if (status != GRB_OPTIMAL && status != GRB_TIME_LIMIT) {

				fprintf(stderr, "Error: it isn't optimal\n");
				goto QUIT;
			}

			error = GRBgetdblattr(modelMIP, GRB_DBL_ATTR_OBJVAL, &objval);
			if (error) goto QUIT;

			error = obtemResultado(modelMIP, &valorVariaveis, m);
			if (error) goto QUIT;

		}

		bestIteration = bestOfBest;

		while (bestIteration) {

			erro = validSoluction(bestIteration->soluction.bins, tamanho, quantidade, bestIteration->soluction.value);

			if (erro > 0) {


				printf("\n\tACO solution not recorded! Found %i errors!\n", erro);

			}
			else {

				saveFile(bestIteration->soluction.bins, classe, quantidade, instancia, bestIteration->soluction.value, bestIteration->soluction.time, semente);

			}

			bestIteration = bestIteration->nextBestSoluction;

		}


		salvaSolucao(&columnSoluction, primeiraColuna, valorVariaveis, m, (int)objval, totalTime);

		excluiItensDuplicados(usedItens, columnSoluction, itens);

		erro = validSoluction(columnSoluction->soluction.bins, tamanho, quantidade, columnSoluction->soluction.value);

		if (erro > 0) {

			printf("\n\tCG solution not recorded! Found %i errors!\n", erro);

		}else{
		
			saveFile(columnSoluction->soluction.bins, classe, quantidade, instancia, columnSoluction->soluction.value, columnSoluction->soluction.time, semente);	
			
		}
		
		printf("\n\tSuccessfully created solutions! \n\tResults recorded in results.csv\n");

		printf("\n\n\t++++++++++++++++++++++++++++++++++++++++++++\n");

		printf("\n\t%d_%d_%d\n", classe, quantidade, instancia);

		printf("\n\n\tType 1 to solve other problem or 0 to exit!\n\t");
		scanf("%d", &opcao);

		for (k = 0; k < quantidade; k++) {

			free(itens[k]);
			free(nextStep[k]);

		}

		free(itens);
		free(ordem);
		free(bestOrder);
		free(usedItens);
		free(nextStep);
		free(PiValue);

	}

	/* Report the result */


QUIT:

	/* Error reporting */



	if (error) {

		printf("\nErro %d\n", error);

		printf("ERROR: %s\n", GRBgeterrormsg(env));
		system("pause");
		exit(1);
	}


	//Free model
	GRBfreemodel(model);

	GRBfreemodel(modelMIP);

	//Free environment
	GRBfreeenv(env);

	return 0;

}
