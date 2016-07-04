// Main file for hybrid compression


#include "huffman.h"
#define SAMPLING_FILE "save_data.txt"
#define BYTES_LEN 16 // for 8 byte data
#define FILE_SIZE 100 // change it inorder to adjust value
#define EXP_FILE "save_exp.txt"
#define MH_FILE "save_mh.txt"
#define ML_FILE "save_ml.txt"

#define SIGN 1
#define EXP 11
#define MH 20
#define ML 32

using namespace std;

// write results in file
void write_file(string file_name,unordered_map<string,string> code)
{
	ofstream ofile;
	ofile.open(file_name);
	auto lt = code.begin();
	for(;lt!=code.end();lt++)
	{
		pair<string,string> p = *lt;
		ofile<<p.first<<" "<<p.second<<endl;
	}
}
// traverse to form act code
void traverse(hufman *root,string stack,unordered_map<string,string>& code)
{
	if(!root) return;

	if(root->is_leaf)
	{
		code.insert(make_pair(root->symbol,stack));
		return;
	}

	stack.push_back('0');
	traverse(root->left,stack,code);
	stack.pop_back();
	stack.push_back('1');
	traverse(root->right,stack,code);
	stack.pop_back();
}
// To split each line -> SKIP in GPUSIM
vector<string> split(string input,char del)
{
	stringstream ss (input);
	string each;
	vector<string> res;
	while(getline(ss,each,del))
	{
		res.push_back(each);
	}
	return res;
}
// to convert hex to decimal
unsigned int get_dec(string raw)
{
	stringstream convert(raw);
	unsigned int ret;
	convert >> std::hex >> ret;
	return ret;
}
// function to get binary string from dec
string get_bin_str(unsigned int val)
{
	string res;
	int count=4; // each 4 bit only , need to change for generic code
	while(count>0)
	{
		res.insert(res.begin(),(val&1)+'0');
		val>>=1;
		count--;
	}
	return res;
}
// function to convert hexadecimal to binary // returns output in form of string(0,1)
vector<string> convert_binary(string raw_data)
{
	vector<string> bin_data;
	for(int i=0;i<raw_data.size();i+=BYTES_LEN)
	{
		string cur = raw_data.substr(i,BYTES_LEN);
		string res;
		for(int j = 0; j<cur.size();j++)
		{
			string tmp;
			tmp.push_back(cur[j]);
			unsigned int val = get_dec(tmp);
			res+=get_bin_str(val);
		}
		bin_data.push_back(res);
	}
	return bin_data;
}
// function to read raw data -> SKIP in GPUSIM
void read_data(string file_name,vector<string>& raw_data)
{
	ifstream ifile;
	ifile.open(file_name);
	string input;
	while(getline(ifile,input))
	{
		vector<string> each_line = split(input,' ');
		raw_data.push_back(each_line[3]);
	}
}
// split and count for each
void split_raw_value(unordered_map<string,int>& exp_track,unordered_map<string,int>& mh_track,
						unordered_map<string,int>& ml_track,vector<string> raw_binary_data)
{
	for(int i=0;i<raw_binary_data.size();i++)
	{
		string cur = raw_binary_data[i];
		string exp = cur.substr(1,EXP);
		string mh = cur.substr(EXP+1,MH);
		string ml = cur.substr(1+EXP+MH);
		
		if(exp_track.find(exp)==exp_track.end())
		exp_track.insert(make_pair(exp,1));
		else
		exp_track[exp]++;
		
		if(mh_track.find(mh)==mh_track.end())
		mh_track.insert(make_pair(mh,1));
		else
		mh_track[mh]++;
		
		if(ml_track.find(ml)==ml_track.end())
		ml_track.insert(make_pair(ml,1));
		else
		ml_track[ml]++;
	}
}
bool myfunc_dec(pair<int,pair<hufman*,string> > a,pair<int,pair<hufman*,string> > b)
{
	return a.first>b.first;
}
bool myfunc_inc(pair<int,pair<hufman*,string> > a,pair<int,pair<hufman*,string> > b)
{
	return a.first<b.first;
}
hufman* construct_tree(vector<pair<int,pair<hufman*,string> > >& track)
{
	if(track.size()==0) return NULL;

	if(track.size()==1 && track[0].second.first!=NULL) return track[0].second.first;

	pair<int,pair<hufman*,string> > a = track[0];
	pair<int,pair<hufman*,string> > b = track[1];
	
	track.erase(track.begin());
	track.erase(track.begin());

	hufman *node = new hufman(false,"");
	node->val=a.first+b.first;
	hufman *one_node = a.second.first;
	hufman *two_node = b.second.first;
	if(!a.second.first)
	{
		one_node = new hufman(true,a.second.second);
		one_node->val=a.first;	
	}
	if(!b.second.first)
	{
		two_node = new hufman(true,b.second.second);
		two_node->val=b.first;
	}
	
	node->left=one_node;
	node->right=two_node;

	pair<int,pair<hufman*,string> > c;
	c.first=a.first+b.first;
	c.second.first=node;

	track.push_back(c);
	sort(track.begin(),track.end(),myfunc_inc);
	
	return construct_tree(track);	
					
}
// void fill data
void fill_data(unordered_map<string,int> src_track,vector<pair<int,pair<hufman*,string> > >& src_tree)
{
	auto it = src_track.begin();
	for(;it!=src_track.end();it++)
	{
		pair<string,int> p = *it;
		pair<hufman*,string> sub;
		sub.first=NULL;
		sub.second.assign(p.first);
		pair<int,pair<hufman*,string> > nor;
		nor.first=p.second;
		nor.second=sub;
		src_tree.push_back(nor);	
	}
	sort(src_tree.begin(),src_tree.end(),myfunc_dec);
	if(src_tree.size()>FILE_SIZE)
	{
		src_tree.resize(FILE_SIZE);
	}
	sort(src_tree.begin(),src_tree.end(),myfunc_inc);
}
// code sampling function
void sampling_data(string file_name)
{

	unordered_map<string,int> exp_track;
	unordered_map<string,int> mh_track;
	unordered_map<string,int> ml_track;
	vector<string> raw_data;
	read_data(file_name,raw_data);
	
	for(int i=0;i<raw_data.size();i++)
	{
		// raw_data is single packet we get in gpusim
		// can eliminate this for loop in act code
		// convert hex value in binary
		vector<string> raw_binary_data = convert_binary(raw_data[i]);
	
		// split value in exponent, mantisa high and mantisa low
		split_raw_value(exp_track,mh_track,ml_track,raw_binary_data);
	}
	
	// convert data to huffman tree
	vector<pair<int,pair<hufman*,string> > > exp_tree;
	vector<pair<int,pair<hufman*,string> > > mh_tree;
	vector<pair<int,pair<hufman*,string> > > ml_tree;
	
	unordered_map<string,string> exp_code;
	unordered_map<string,string> mh_code;
	unordered_map<string,string> ml_code;
	
	fill_data(exp_track,exp_tree);
	hufman *exp_hufman = construct_tree(exp_tree);
	traverse(exp_hufman,"",exp_code);
	write_file(EXP_FILE,exp_code);
	
	fill_data(mh_track,mh_tree);
	hufman *mh_hufman = construct_tree(mh_tree);
	traverse(mh_hufman,"",mh_code);
	write_file(MH_FILE,mh_code);
	
	fill_data(ml_track,ml_tree);
	hufman *ml_hufman = construct_tree(ml_tree);
	traverse(ml_hufman,"",ml_code);
	write_file(ML_FILE,ml_code);
}
void init_tracker(string file_name,unordered_map<string,string>& src_tracker)
{
	ifstream ifile;
	ifile.open(file_name);
	string input;
	while(getline(ifile,input))
	{
		vector<string> split_data = split(input,' ');
		src_tracker.insert(make_pair(split_data[0],split_data[1]));
	}
}
string do_encode(vector<string> input,float& orig_bytes,float& encode_bytes,unordered_map<string,string> exp_tracker,
								unordered_map<string,string> mh_tracker,unordered_map<string,string> ml_tracker)
{
	string output;
	// input refers to single packet, i.e 128 bytes packet split into 8 bytes words
	
	for(int i=0;i<input.size();i++)
	{
		string cur = input[i];
		string exp = cur.substr(1,EXP);
		string mh = cur.substr(EXP+1,MH);
		string ml = cur.substr(1+EXP+MH);
		orig_bytes+=(cur.size()/8.0);
		string exp_code = exp_tracker.find(exp)!=exp_tracker.end() ? exp_tracker[exp] : exp; 
		string mh_code = mh_tracker.find(mh)!=mh_tracker.end() ? mh_tracker[exp] : mh;
		string ml_code = ml_tracker.find(ml)!=ml_tracker.end() ? ml_tracker[exp] : ml;
		
		output = output + cur.substr(0,1) + exp_code + mh_code + ml_code; 
	}
	encode_bytes+=(output.size()/8.0);
	return output;
}
void encode_data(string file_name)
{
	// this part can be skipped in gpgpusim
	vector<string> raw_data;
	read_data(file_name,raw_data);
	
	// init trackers;
	unordered_map<string,string> exp_tracker;
	unordered_map<string,string> mh_tracker;
	unordered_map<string,string> ml_tracker;
	
	init_tracker(EXP_FILE,exp_tracker);
	init_tracker(MH_FILE,mh_tracker);
	init_tracker(ML_FILE,ml_tracker);
	
	float orig_bytes=0,encode_bytes=0;
	for(int i=0;i<raw_data.size();i++)
	{
		vector<string> raw_binary_data = convert_binary(raw_data[i]);
		// single packet is divided into 16 bytes raw_binary data
		string output = do_encode(raw_binary_data,orig_bytes,encode_bytes,exp_tracker,mh_tracker,ml_tracker);
	}
	cout<<"Original Bytes : "<<orig_bytes<<" Encode Bytes : "<<encode_bytes<<endl;
}
int main(int argc,char *argv[])
{
	if(argc<2)
	{
		cout<<" Please enter file name , Sampling (0) || Encoding (1)" << endl;
		return 0;
	}
	
	if(stoi(argv[2])==0)
	{
		sampling_data(argv[1]);
	}
	else
	{
		encode_data(argv[1]);
	}
	return 0;
}