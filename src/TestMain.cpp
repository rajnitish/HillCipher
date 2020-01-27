#include "HillCipher.h"

int main() {


	char ch;
	while(1)
	{
		cout<<"\n *************** HILL CIPHER ******************"<<endl;
		cout<<"Type E for Encrypt or D for Decrypt C for CryptAnalysis. Press Esc for exit"<<endl;

		ch = getch();

		string ip_file,key_file,op_file;
		HillCipher hc_e,hc_d,hc_c;
		switch(ch)
		{
		case 'E':
		{
			cout<<"Enter Key File Name"<<endl;
			cin>>key_file;

			if(!hc_e.ReadKeyMatrix(key_file)) 	break;

			cout<<" Enter Plain Text File Name"<<endl;
			cin>>ip_file;
			if(!hc_e.ReadPlainText(ip_file)) 	break;

			hc_e.performEncryption();

			break;
		}
		case 'D':
		{
			cout<<"Enter Key File Name"<<endl;
			cin>>key_file;

			if(!hc_d.ReadKeyMatrix(key_file)) 	break;

			hc_d.PerformInverseKeyMatrix();

			hc_d.pad = hc_e.pad;
			cout<<" Enter Cipher Text File Name"<<endl;
			cin>>ip_file;
			if(!hc_d.ReadCipherText(ip_file)) 	break;

			hc_d.performDecryption();


			break;
		}
		case 'C':
			cout<<" Enter Cipher Text File Name"<<endl;
			cin>>ip_file;
			if(!hc_c.ReadCipherText(ip_file)) 	break;
			hc_c.performCryptAnalysis();
			break;
		case 27:
			cout<<" Thanks!!! "<<endl;
			sleep(4);
			exit(0);
		default:
		{
			cout<<" Try Again!!! Wrong input"<<endl;
			break;
		}

		}
	}

}
