#ifndef _Individual_
#define _Individual_

class Individual {

private:

protected:
    int m_N; 	//Number of genes (cities)
    int* m_s; 	//Sequence of genes 
   // double fitness;

public:
    // constructors
    Individual();  //sets the sequence as [1 2 3 ....32]
    Individual(int N, int* x); //sets a specified sequence
    // destructor
    ~Individual();
    // methods
    int* GetSequence(); //returns the sequence of genes
    int GetGene(int i); //returns the value of the gene in posiion i. First gene is at i = 0
    void SetSequence(int* x);
    void SetGene(int i, int x);
};

#endif 
