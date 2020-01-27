#ifndef __HILLCIPHER_CPP__
#define __HILLCIPHER_CPP__

#include "HillCipher.h"

twoDArray HillCipher::PerformInverse(int KEY[N][N],int size,bool &isInvertible)
{
	twoDArray obj;

	int det;
	if (size == 1)	det =  KEY[0][0];
	else		det = fmod((ComputeDeterminant(KEY, size)+26),26);

	if (det == 0)
	{
		//cout << "Not Invertible.Singular Matrix"<<endl;
		isInvertible = false;
		return obj;
	}
	if( __gcd(det, 26) != 1)
	{
		//cout << "Not Invertible. Determinant is not Co-Prime"<<endl;
		isInvertible = false;
		return obj;

	}

	int AdjointMatrix[N][N];
	ComputeAdjoint(KEY, AdjointMatrix);

	int mulInvDet = 1;
	det = det%26;
	for (; mulInvDet<26; mulInvDet++)
		if ((det*mulInvDet) % 26 == 1)
			break;

	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			obj.A[i][j] = fmod((AdjointMatrix[i][j])*mulInvDet, 26);
			obj.A[i][j] = fmod((obj.A[i][j]+26),26);

		}
	}

	isInvertible = true;

	return obj;
}

void HillCipher::printMatrix(int M[N][N],int size)
{
	for(unsigned char i = 0 ;i<size;i++ )
	{
		for(unsigned char j = 0 ;j<size;j++ )
		{
			cout<<M[i][j]<<" ";

		}
		cout<<endl;
	}
}

twoDArray HillCipher::matrixMul(int A[N][N],int B[N][N],int size)
{
	twoDArray obj;

	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			obj.A[i][j]=0;
			for(int k=0;k<size;k++)
			{
				obj.A[i][j] += fmod(A[i][k],26) * fmod(B[k][j],26);
			}
			obj.A[i][j] = fmod(obj.A[i][j], 26);
		}
	}

	return obj;

}
float ComputeIndexOfCoincidence(vector<unsigned char> file)
{
	int freqEachChar[26] = {0};
	float IC = 0.0;
	int NCHAR = file.size();
	for(int i = 0;i <NCHAR;i++)
	{
		freqEachChar[file[i]-97]++;
	}

	for(int i = 0 ; i<26;i++)
		IC += (freqEachChar[i]) *( (freqEachChar[i]) -1 );


	IC = 26*IC/(NCHAR*(NCHAR-1));

	return IC;
}
void HillCipher::performCryptAnalysis()
{
	ParseCipherForDigram();

	int matForInv[N][N];

	M_size = 2;
	int KEY[N][N];
	int K = 1;
	for(int i = 0 ; i<Digraph.size()-1;i++)
		for(int j = i+1 ;j<Digraph.size();j++)
		{
			string str = Digraph[i] + Digraph[j];

			for(unsigned char I = 0 ;I<M_size;I++ )
			{
				for(unsigned char J = 0 ;J<M_size;J++ )
				{
					matForInv[I][J] = str[M_size*I + J];
					matForInv[I][J] = matForInv[I][J] - 97;
					//cout<<matForInv[i][j]<<" ";
				}
				//cout<<endl;
			}
			bool isInvert;
			twoDArray INVTD_MTRX =	PerformInverse(matForInv,M_size,isInvert);

			if(!isInvert)continue;

			int FREQ_GRAM[N][N];
			for(int C = 0 ; C <Cipher_Digram.size()-1;C++)
				for(int C1 = C+1 ; C1 <Cipher_Digram.size();C1++)
				{
					string str = Cipher_Digram[C].dia + Cipher_Digram[C1].dia;
					if(Cipher_Digram[C].freq<=1 || Cipher_Digram[C1].freq <=1) continue;
					for(unsigned char ii = 0 ;ii<M_size;ii++ )
					{
						for(unsigned char jj = 0 ;jj<M_size;jj++ )
						{
							FREQ_GRAM[ii][jj] = str[M_size*ii + jj];
							FREQ_GRAM[ii][jj] = FREQ_GRAM[ii][jj] - 97;
							//cout<<FREQ_GRAM[ii][jj]<<" ";

						}
						//	cout<<endl;
					}

					twoDArray GEN_KEY;
					GEN_KEY = matrixMul(FREQ_GRAM,INVTD_MTRX.A,M_size);

					//if (GEN_KEY.A[0][0] == 22 && GEN_KEY.A[1][1] == 6 && GEN_KEY.A[0][1] == 3 && GEN_KEY.A[1][0] == 9)
					{

						PlainTextVector.clear();
						int i1, j1, k1;
						for(j1 = 0; j1 < CipherTextVector.size()/M_size; j1++) // traversing Cipher text vector and treat as row major vector
						{
							for(i1 = 0; i1 < M_size; i1++) // Column no of Key Matrix
							{
								int tm = 0;
								for(k1 = 0; k1 < M_size; k1++) // Row no of Key Matrix
								{
									bool flag = false;
									twoDArray InvGEN_KEY = PerformInverse(GEN_KEY.A,M_size,flag);
									if(!flag) break;

									tm = tm + fmod(InvGEN_KEY.A[i1][k1],26) * fmod(CipherTextVector[j1*M_size + k1] -97,26);
								}

								tm = fmod(tm, 26);
								PlainTextVector.push_back((unsigned char)(tm+97));
							}
						}

						float IC = ComputeIndexOfCoincidence(PlainTextVector);

						if(IC>1.73&& IC<26){
							cout<<" IC = "<<IC<<endl;
							out_file.open("Plain_PostCryAnalysis_"+to_string(K)+".txt", ios::out | std::ios::binary );
							K++;
							if (out_file.fail()) {
								cout << "Couldn't open the file! :: FILE NOT EXIST \n" << endl;
								return;
							}

							for(int m = 0; m < PlainTextVector.size(); m++)
								out_file<<PlainTextVector[m];

							out_file.close();
						}


					}

				}

		}

}



bool sortbyfreq(const struct_digram &a,const struct_digram &b)
{
	return (a.freq>b.freq);
}

void HillCipher::ParseCipherForDigram()
{
	//cout<<" CipherTextVector.size()-1 "<<CipherTextVector.size()-1<<endl;
	for(int i = 0; i<CipherTextVector.size()-1;i++)
	{
		string str1(1,CipherTextVector[i]);
		string str2(1,CipherTextVector[i+1]);

		updateDigram(str1+str2);
		//cout<<" "<<str1+str2;
	}
	//cout<<endl;

	sort(Cipher_Digram.begin(), Cipher_Digram.end(), sortbyfreq);

	//cout<<"After sorting --------->"<<endl;
	for(int i = 0 ; i <Cipher_Digram.size();i++)
	{
		//Cipher_Digram[i].print();
	}
}
void HillCipher::updateDigram(string ch)
{
	for(int i =0 ; i<Cipher_Digram.size();i++)
	{
		if(ch == Cipher_Digram[i].dia)
		{
			Cipher_Digram[i].freq++;
			return;
		}
	}

	struct_digram obj;
	obj.freq = 1;
	obj.dia = ch;
	Cipher_Digram.push_back(obj);

}

bool HillCipher::ReadPlainText(string ip_file)
{
	std::string str = ip_file;
	in_file.open(str, ios::in | std::ios::binary );

	if (!in_file) {
		cout << "Couldn't open the file! :: FILE NOT EXIST \n" << endl;
		return false;
	}

	unsigned char P;
	unsigned int elem_count = 0;
	while (!in_file.eof())
	{
		P = in_file.get();
		if(P<=122 && P>= 97)
		{
			PlainTextVector.push_back(P);
			elem_count++;
		}
	}

	in_file.close();
	if(0 == elem_count)
	{
		cout<<" Plain text File is not compatible for Hill Cipher :: No Element present \n"<<endl;
		return false;
	}
	else if(elem_count%M_size != 0)
	{
		cout<<" Plain text File is not compatible for Hill Cipher :: #Elements is not multiple of KeyMatrix Order  \n"<<endl;
		cout<<" Would you like to padding with random characters?? Press Y else Exit ... Decrypter must know this info to avoid these bytes"<<endl;
		char k;
		cin>>k;
		if(k == 'Y')
		{
			pad = (M_size - elem_count%M_size);
			for(int i = 0 ;i < (M_size - elem_count%M_size) ; i++ )
			{
				PlainTextVector.push_back('a');
			}
		}
		else
			return false;
	}
	/*
	cout<<" CHECKING ELEMENT \n\n"<<PlainTextVector.size()<<endl;

	for(unsigned char i = 0 ;i<PlainTextVector.size();i++ )
	{
		cout<<(int)PlainTextVector[i]<<" " ;

	}
	 */

	return true;


}
bool HillCipher::ReadKeyMatrix(string key_file)
{

	std::string str = key_file;
	in_file.open(str, ios::in | std::ios::binary );

	if (in_file.fail()) {
		cout << "Couldn't open the file! :: FILE NOT EXIST \n" << endl;
		return false;
	}

	unsigned int elem_count = 0;
	int a;
	vector<int> elem;
	while (in_file >> a)
	{
		elem_count++;
		elem.push_back(a);
	}
	in_file.close();

	if(0 == elem_count)
	{
		cout<<" Key File is not compatible for Hill Cipher :: No Element present KEY MATRIX  \n"<<endl;
		return false;
	}
	else if(std::floor(sqrt(elem_count)) != sqrt(elem_count))
	{
		cout<<" Key File is not compatible for Hill Cipher :: No Perfect Square KEY MATRIX  \n"<<endl;
		return false;
	}

	M_size = sqrt(elem_count);

	for(unsigned char i = 0 ;i<M_size;i++ )
	{
		vector<int> tmp;
		for(unsigned char j = 0 ;j<M_size;j++ )
		{
			tmp.push_back(elem[M_size*i + j]);
		}

		KeyMatrix.push_back(tmp);

	}
	/*
	cout<<" CHECKING ELEMENT \n\n"<<endl;
	for(unsigned char i = 0 ;i<M_size;i++ )
		{
			for(unsigned char j = 0 ;j<M_size;j++ )
			{
				cout<< KeyMatrix[i][j]<<" ";
			}
			cout<<endl;
		}
	 */

	return true;


}

void HillCipher::performEncryption()
{

	int i, j, k;


	for(j = 0; j < PlainTextVector.size()/M_size; j++) // traversing plain text vector and treat as row major vector
	{
		for(i = 0; i < M_size; i++) // Column no of Key Matrix
		{
			int tm = 0;
			for(k = 0; k < M_size; k++) // Row no of Key Matrix
			{
				tm = tm + fmod(KeyMatrix[i][k],26) * fmod(PlainTextVector[j*M_size + k] -97,26); // Considering Sir has given key coloumn wise
				//tm = tm + fmod(KeyMatrix[k][i],26) * fmod(PlainTextVector[j*M_size + k] -97,26);
			}

			tm = fmod(tm, 26);
			CipherTextVector.push_back((unsigned char)(tm+97));
		}
	}

	out_file.open("Cipher.txt", ios::out | std::ios::binary );

	if (out_file.fail()) {
		cout << "Couldn't open the file! :: FILE NOT EXIST \n" << endl;
		return;
	}

	for(i = 0; i < CipherTextVector.size(); i++)
		out_file<<CipherTextVector[i];

	out_file.close();
}

void HillCipher::ComputeCofactor(int A[N][N], int temp[N][N], int ROW, int COL, int n)
{
	int i = 0, j = 0;

	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++)
		{
			if (r != ROW && c != COL)
			{
				temp[i][j++] = A[r][c];
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}

}

int HillCipher::ComputeDeterminant(int A[N][N], int n)
{
	int Det = 0;
	if (n == 1)		return A[0][0];
	int coFac[N][N];
	int sign = 1;  // To store sign multiplier

	for (int i = 0; i < n; i++)
	{
		ComputeCofactor(A, coFac, 0, i, n);
		Det += sign * A[0][i] * ComputeDeterminant(coFac, n - 1);
		sign = -sign;
	}

	return fmod((Det+26),26);
}

void HillCipher::ComputeAdjoint(int A[N][N],int adj[N][N])
{
	if (N == 1)
	{
		adj[0][0] = 1;
		return;
	}

	int sign = 1, temp[N][N];

	for (int i=0; i<M_size; i++)
	{
		for (int j=0; j<M_size; j++)
		{
			ComputeCofactor(A, temp, i, j, M_size);
			sign = ((i+j)%2==0)? 1: -1;
			adj[j][i] = (sign)*(ComputeDeterminant(temp, M_size-1));
		}
	}
}

void HillCipher::PerformInverseKeyMatrix()
{
	int KEY[N][N],inverse[N][N];

	for(int i = 0; i< M_size;i++)
		for(int j = 0 ;j <M_size ;j++)
			KEY[i][j] = KeyMatrix[i][j];

	int det;
	if (M_size == 1)	det =  KEY[0][0];
	else		det = fmod((ComputeDeterminant(KEY, M_size)+26),26);

	if (det == 0)
	{
		cout << "Not Invertible."<<endl;
		return;
	}

	int AdjointMatrix[N][N];
	ComputeAdjoint(KEY, AdjointMatrix);

	int mulInvDet = 1;
	det = det%26;
	for (; mulInvDet<26; mulInvDet++)
		if ((det*mulInvDet) % 26 == 1)
			break;

	for (int i=0; i<M_size; i++)
	{
		vector<int> tmp;
		for (int j=0; j<M_size; j++)
		{
			inverse[i][j] = fmod((AdjointMatrix[i][j])*mulInvDet, 26);
			inverse[i][j] = fmod((inverse[i][j]+26),26);
			tmp.push_back(inverse[i][j]);
		}
		InvKeyMatrix.push_back(tmp);
	}

}
void HillCipher::performDecryption()
{
	int i, j, k;
	for(j = 0; j < CipherTextVector.size()/M_size; j++) // traversing Cipher text vector and treat as row major vector
	{
		for(i = 0; i < M_size; i++) // Column no of Key Matrix
		{
			int tm = 0;
			for(k = 0; k < M_size; k++) // Row no of Key Matrix
			{
				tm = tm + fmod(InvKeyMatrix[i][k],26) * fmod(CipherTextVector[j*M_size + k] -97,26); // Considering Sir has given key coloumn wise
				//tm = tm + fmod(InvKeyMatrix[k][i],26) * fmod(CipherTextVector[j*M_size + k] -97,26);
			}

			tm = fmod(tm, 26);
			PlainTextVector.push_back((unsigned char)(tm+97));
		}
	}

	out_file.open("Plain_PostDecrypt.txt", ios::out | std::ios::binary );

	if (out_file.fail()) {
		cout << "Couldn't open the file! :: FILE NOT EXIST \n" << endl;
		return;
	}

	for(i = 0; i < PlainTextVector.size(); i++)
		out_file<<PlainTextVector[i];

	out_file.close();

}

bool HillCipher::ReadCipherText(string ip_file)
{
	std::string str = ip_file;
	in_file.open(str, ios::in | std::ios::binary );

	if (!in_file) {
		cout << "Couldn't open the file! :: FILE NOT EXIST \n" << endl;
		return false;
	}

	unsigned char P;
	unsigned int elem_count = 0;
	while (!in_file.eof())
	{
		P = in_file.get();
		if(P<=122 && P>= 97)
		{
			CipherTextVector.push_back(P);
			elem_count++;
		}
	}
	in_file.close();
	if(0 == elem_count)
	{
		cout<<" Cipher text File is not compatible for Hill Cipher :: No Element present \n"<<endl;
		return false;
	}

	for(int i = 0; i<pad;i++)
		CipherTextVector.pop_back(); //Since padding was done while reading plain text

	return true;

}

#endif
