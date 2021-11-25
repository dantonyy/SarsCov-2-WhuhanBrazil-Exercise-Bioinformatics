# coding: utf-8

import re
import sys
import Bio # Biblioteca Biophyton, necessária para tratamento da sequência genética
from Bio.Seq import Seq #Funçao .Seq responsável por obter o sequenciamento genético, reverso complementar, RNA e Proteina
import difflib # Biblioteca difflib, necessária para realização da função de diferenciação entre strings

#Abre o arquivo do gene S de Whuhan como dna_whu
with open("gene_S_Wuhan.txt", "r") as arquivo:
    dna_whu = arquivo.read()

#Abre o arquivo do gene S do brazil como dna_br
with open("gene_S_Brazil.txt", "r") as arquivo_br:
    dna_br = arquivo_br.read()

# Define DNA e RNA para posteriormente verificar proporções dos mesmos
dna = 'ATGC'
rna = 'AUGC'

#Menu
print('''[ 1 ] Verificar as proporções de A T G C na sequência de Whuhan
[ 2 ] Verificar a sequência proteica 
[ 3 ] Inserir curta sequência de nucleotídeos, e comparar se está presente no gene S de Whuhan, apontando indice do começo da sequencia
[ 4 ] Mostrar principais alterações e index das alterações dos genes de Whuhan e Brazil
[ 5 ] Sair do programa''')
opcao = 0
while opcao != 5:
    opcao = int(input('Opção desejada: '))
    if opcao == 1:
        # Responde a Questão 1 : Verificar a proporção de C e G (GC content) na sequência

        #Linha 44~47: Percorre toda cadeia de caracteres da variavel dna_whu e faz a contagem de caracteres armazenados na variável dna.
        dna_count = {}.fromkeys(dna, 0)
        for w in dna_whu:
            if w in dna_count:
                dna_count[w] += 1
        #Linha 49, 50: Imprime proporções de A T G C na sequência de Whuhan
        print("1. Verificando proporções de A T G C na sequência:")
        print(dna_count)

    if opcao == 2:
        # Responde a Questão 2 : Converter a sequência de nucleotídeos na respectiva sequência proteica

        # Linha 56~60: Efetua a transcrição do DNA em RNA, fazendo a substituição por seus pares de bases complementares
        transc = dna_whu.replace("A", "U")
        transc2 = transc.replace("T", "A")
        transc3 = transc2.replace("G", "H")
        transc4 = transc3.replace("C", "G")
        transc5 = transc4.replace("H", "C")

        # Linha 63~66: Percorre toda cadeia de caracteres da variavel transc5 (RNA) e faz a contagem de caracteres armazenados na variável rna.
        transc_count = {}.fromkeys(rna, 0)
        for w in transc5:
            if w in transc_count:
                transc_count[w] += 1

        # Linha 69~71: Imprime proporções de A U G C no RNAm
        print("\n2. Converter a sequência de nucleotídeos na respectiva sequência proteica")
        print("\nProporções de A U G C Após a transcrição: ")
        print(transc_count, "\n")

        whu_seq1 = Seq(dna_whu) # Utilizando a biblioteca Biophyton, faz o sequenciamento do DNA
        whu_protein = whu_seq1.translate() # Transforma o DNA obtido anteriormente em proteina

        print("A sequência da proteína é", whu_protein) # Imprime sequência proteica

    if opcao == 3:
        # Responde a Questão 3 : Verificar se uma outra sequência curta de nucleotídeos, passada pelo usuário do seu sistema, está presente no gene S e, se positivo, indicar o índice do início da sequência.

        print("\n3. Verificar se uma outra sequência curta de nucleotídeos, passada pelo usuário do seu sistema, está presente no gene S e, se positivo, indicar o índice do início da sequência.")

        usuario_nuc = input("\nInsira sequência curta de nucleotídeos:")

        print("\nÍndice do primeiro item da sequencia: \n", dna_whu.find(usuario_nuc))

    if opcao == 4:
        # Respnode a Questão 4 : Comparar a sequência de nucleotídeos com a sequência de Wuhan;
        # Responde a Questão 5 : Comparar a sequenciado proteica com a sequência de Wuhan;
        # Responde a Questão 6 : Mostrar as principais alterações (qual nucleotídeo / aminoácido foi substituído e qual a posição / index da substituição)
        br_seq = Seq(dna_br)
        br_protein = br_seq.translate()

        # Linhas 95~97: Transforma a sequência corrida em trincas de base para melhor comparação das alterações nos nucleotídeos e aminoácidos
        N = 3
        list_whu = [dna_whu[n:n + N] for n in range(0, len(dna_whu), N)]
        list_br = [dna_br[n:n + N] for n in range(0, len(dna_br), N)]

        # Utilizando a biblioteca difflib, imprime toda sequência de trincas de bases e aponta as principais alterações nos nucleotideos e aminoácidos
        diff = difflib.unified_diff(list_whu, list_br)
        print('\n'.join(list(diff)))