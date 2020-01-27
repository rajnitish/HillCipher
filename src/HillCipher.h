#ifndef __HILLCIPHER_H__
#define __HILLCIPHER_H__

#include<iostream>
#include<math.h>
#include <vector>
#include <fstream>
#include <string>
#include <conio.h>
#include <unistd.h>
#include<bits/stdc++.h>

using namespace std;

#define N 10

struct struct_digram
{
	string dia;
	int freq;

	void print()
	{
		cout<<" GRAM = "<<dia<<" FREQ = "<<freq<<endl;
	}
};

struct twoDArray
{
	int A[N][N];
};
class HillCipher
{
	unsigned char M_size;
	vector< vector<int> > KeyMatrix;
	vector< vector<int> > InvKeyMatrix;
	vector<unsigned char> PlainTextVector;
	vector<unsigned char> CipherTextVector;
	fstream in_file, out_file;
	vector<string> Digraph = {"th", "er", "on", "an" ,"re", "he", "in", "ed", "nd", "ha", "at", "en",
			                     "es", "of", "or", "nt", "ea",  "ti", "to", "it", "st", "io", "le", "is","ou", "ar", "as", "de", "rt", "ve"};

	vector<struct_digram> Cipher_Digram;
public:

	int pad = 0;
	HillCipher(){M_size = 0;}
	bool ReadKeyMatrix(string);
	bool ReadPlainText(string);
	void performEncryption();
	void PerformInverseKeyMatrix();
	void performDecryption();
	int ComputeDeterminant(int A[N][N], int);
	void ComputeAdjoint(int A[N][N],int adj[N][N]);
	void ComputeCofactor(int A[N][N], int temp[N][N], int, int, int);
	bool ReadCipherText(string);
	void updateDigram(string);
	void ParseCipherForDigram();
	void performCryptAnalysis();
	twoDArray PerformInverse(int KEY[N][N],int,bool&);
	void printMatrix(int M[N][N],int size);
	twoDArray matrixMul(int A[N][N],int B[N][N],int size);

};

#endif
