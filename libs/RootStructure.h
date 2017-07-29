#ifndef REPRESENTATION_LIEALGEBRA_H
#define REPRESENTATION_LIEALGEBRA_H


#include<Eigen/Core>
#include<vector>
#include<boost/rational.hpp>
#include<boost/type_index.hpp>
#include<iostream>
#include<algorithm>
#include<memory>
#include"LoggerLibrary.h"

// Tells Eigen how to handle type
namespace Eigen {
	template<> struct NumTraits<boost::rational<int>> :GenericNumTraits<boost::rational<int>> {
		enum {
			IsInteger = 1,
		};
	};
}

namespace Representation {
	using namespace boost;
	using namespace Eigen;

	/*
	* Allows for more compact handling of roots and weights
	*/
	struct weight {
		typedef Matrix <rational<int>, Dynamic, 1> VectorXr;

		VectorXr alpha;
		VectorXr omega;
		VectorXr ortho;

		bool is_empty();

		bool operator==(const weight& w);
		weight& operator=(const weight& w);
		weight operator+(const weight& rhs);
		weight operator*(int i);


		weight() {}
		weight(VectorXr omega);
	};

	weight operator*(int i, const weight& rhs);
	weight operator*(const weight& rhs, int i);



	/*
	*           Helper functions
	*/

	template<typename T>
	bool is_pos(Matrix<T, Dynamic, 1> V)
	{
		FILE_LOG(Diagnositics::logDEBUG4) << "Is Positive called with type " << boost::typeindex::type_id<T>().pretty_name();
		for (int i = 0; i < V.size(); i++)
		{
			if (V(i) < 0)
				return false;
		}
		return true;
	}

	template<typename T>
	bool has_zero(Matrix<T, Dynamic, 1> V)
	{
		FILE_LOG(Diagnositics::logDEBUG4) << "has_zero called with type " << boost::typeindex::type_id<T>().pretty_name();

		for (int k = 0; k < V.size(); k++)
		{
			if (V(k) == 0)
				return true;
		}

		return false;
	}

	MatrixXi int_cast_Matrix(Matrix<rational<int>, Dynamic, Dynamic> X);

	Matrix<rational<int>, Dynamic, Dynamic> Identity(int size);

	template<typename T>
	Matrix<T, Dynamic, Dynamic> ReflectionMatrix(Matrix<T, Dynamic, 1> V)
	{
		FILE_LOG(Diagnositics::logDEBUG4)<< "Reflection Matrix with type " <<boost::typeindex::type_id<T>().pretty_name();
		Matrix<T, Dynamic, Dynamic> M;
		M.resize(V.size(), V.size());
		Matrix<T, Dynamic, Dynamic> I = Identity(V.size());
		for (int i = 0; i < V.size(); i++)
		{
			for (int j = 0; j < V.size(); j++)
			{
				M(i, j) = I(i, j) - 2 * V(i)*V(j) / (V.dot(V));
			}
		}
		return M;

	}

	VectorXi int_cast_Vector(Matrix<rational<int>, Dynamic, 1> X);

	template<typename T>
	T RationalDeterminant(Matrix<T, Dynamic, Dynamic> M)
	{
		FILE_LOG(Diagnositics::logDEBUG4)<< "Rational Determinant with type " <<boost::typeindex::type_id<T>().pretty_name();
		int n = M.rows(); // Size of NxN
		T det = 0; //init det

		Matrix<T, Dynamic, Dynamic> m;

		assert(n >= 1 && "An undefined matrix has been passed to RationalDeterminant.");
		if (n == 1) det = M(0, 0);
		else if (n == 2) det = M(1, 1) * M(0, 0) - M(1, 0)*M(0, 1);
		else
		{
			det = 0;

			for (int j1 = 0; j1 < n; j1++)
			{
				m.resize(n - 1, n - 1);
				for (int i = 1; i < n; i++)
				{
					int j2 = 0;
					for (int j = 0; j < n; j++)
					{
						if (j == j1) continue;
						m(i - 1, j2) = M(i, j);
						j2++;
					}
				}
				int pmone = 0;
				if ((2 + j1) % 2 == 0) pmone = 1;
				else pmone = -1;
				det += pmone * M(0, j1) * RationalDeterminant<T>(m);
			}
		}
		return det;
	}

	template<typename T>
	Matrix<T, Dynamic, Dynamic> CoFactor(Matrix<T, Dynamic, Dynamic> M)
	{
		FILE_LOG(Diagnositics::logDEBUG4)<< "CoFactor with type " <<boost::typeindex::type_id<T>().pretty_name();

		int n = M.rows(); // Size of NxN
		Matrix<T, Dynamic, Dynamic> m, c;

		c.resize(n - 1, n - 1);
		m.resize(n, n);

		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
				int i1 = 0;
				for (int ii = 0; ii < n; ii++) {
					if (ii == i) continue;

					int j1 = 0;
					for (int jj = 0; jj < n; jj++)
					{
						if (jj == j) continue;
						c(i1, j1) = M(ii, jj);
						j1++;
					}
					i1++;
				}

				T det = RationalDeterminant<T>(c);
				int pmone;
				if ((2 + i + j) % 2 == 0) pmone = 1;
				else pmone = -1;

				m(i, j) = pmone * det;
			}
		}
		return m;
	}

	template<typename T>
	Matrix<T, Dynamic, Dynamic> RationalInverse(Matrix<T, Dynamic, Dynamic> M)
	{
		FILE_LOG(Diagnositics::logDEBUG4)<< "RationalInverse with type " <<boost::typeindex::type_id<T>().pretty_name();

		Matrix<T, Dynamic, Dynamic> co = CoFactor<T>(M);
		T det = RationalDeterminant<T>(M);
		return co.transpose() / det;
	}

	template<typename T>
	Matrix<T, Dynamic, Dynamic> PseudoInverse(Matrix<T, Dynamic, Dynamic> M)
	{
		FILE_LOG(Diagnositics::logDEBUG4)<< "PseudoInverse with type " <<boost::typeindex::type_id<T>().pretty_name();

		Eigen::Matrix<T, Dynamic, Dynamic> Sq_M = M*M.transpose();
		return M.transpose()*(RationalInverse<T>(Sq_M));
	}

	template<typename T>
	T master_formula(Matrix<T, Dynamic, 1> U, Matrix<T, Dynamic, 1> V)
	{
		FILE_LOG(Diagnositics::logDEBUG4) << "master_formula with type " <<boost::typeindex::type_id<T>().pretty_name();
		T numerator = 2 * U.dot(V);
		T denominator = U.dot(U);
		return numerator / denominator;
	}

	/*
	*           Classical Groups
	*/

	enum class GroupType {
		A,
		B,
		C,
		D,
		G,
		F
	};

	/*
	*           Base class that handles all Lie algebra
	*/
	class GroupBase {
	protected:
		virtual void createOrtho() = 0;
		virtual void createMatrices() = 0;
		virtual void createAllBases() = 0;

	public:
		typedef Matrix <rational<int>, Dynamic, Dynamic> MatrixXr;
		typedef Matrix <rational<int>, Dynamic, 1> VectorXr;

		virtual VectorXr to_alpha(weight V) = 0;
		virtual VectorXr to_omega(weight V) = 0;
		virtual VectorXr to_ortho(weight V) = 0;


		virtual std::vector<weight> weightTower() = 0;
		virtual std::vector<weight> weightTower(weight w) = 0;
		virtual std::vector<weight> tensorProductDecomp(weight w1, weight w2) = 0;

		virtual size_t get_Rank() = 0;
		virtual MatrixXr get_Cartan() const = 0;
		virtual MatrixXr get_Omega() const = 0;
		virtual MatrixXr get_CoCartan() const = 0;
		virtual MatrixXr get_QuadraticForm() const = 0;

		virtual std::vector<weight> get_simple() const = 0;
		virtual std::vector<weight> get_positiver() const = 0;
		virtual std::vector<weight> get_fweight() const = 0;

		virtual int dim(weight w) = 0;
		virtual int k_lvl(weight w) = 0;

		GroupBase() {}
		GroupBase(const GroupBase&) = default;
		GroupBase(GroupBase &&) = default;
		virtual ~GroupBase(){}
		virtual GroupBase* clone() const = 0;
		GroupBase& operator=(const GroupBase&) = default;
		GroupBase& operator=(GroupBase&&) = default;
	};

	template<GroupType T>
	class LieBase: public GroupBase {
		typedef Matrix <rational<int>, Dynamic, Dynamic> MatrixXr;
		typedef Matrix <rational<int>, Dynamic, 1> VectorXr;

		GroupType Group;

		VectorXr rho;

		std::vector<weight> weylOrbit(weight head);
		std::vector<weight> weylOrbit(weight head, std::vector<int> stabilizer);
		std::vector<weight> dominantWeights(std::vector<weight> representation);

		void createOrtho();
		void createMatrices();
		void createAllBases();

		weight chamberRotate(weight w, int& counter);

		int freudenthalsRecursion(weight current, std::vector<weight> domWeights, std::vector<std::pair<int, weight>> stabilizedOrbits);

		std::vector<std::pair<int, weight>> multiplicity(weight highest);

		size_t Rank;
		MatrixXr Cartan;
		MatrixXr Omega;
		MatrixXr CoCartan;
		MatrixXr QuadraticForm;

		std::vector<weight> simple, positiver, fweight;;

	public:

		VectorXr to_alpha(weight V);
		VectorXr to_omega(weight V);
		VectorXr to_ortho(weight V);

		LieBase (){}
		LieBase(const size_t Rank);
 		virtual ~LieBase() {};
		virtual LieBase * clone() const{ return new LieBase(*this); }
		LieBase(const LieBase&) = default;
		LieBase(LieBase&&) = default;
		LieBase& operator=(LieBase&&) = default;


		std::vector<weight> get_simple() const;
		std::vector<weight> get_positiver() const;
		std::vector<weight> get_fweight() const;
		std::vector<weight> weightTower(), weightTower(weight w);
		std::vector<weight> tensorProductDecomp(weight w1, weight w2);

		size_t get_Rank();
		MatrixXr get_Cartan() const;
		MatrixXr get_Omega() const;
		MatrixXr get_CoCartan() const;
		MatrixXr get_QuadraticForm() const;

		int dim(weight w);
		int k_lvl(weight w);

	};

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::weylOrbit(weight head)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "weylOrbit called";

		std::vector<weight> master_list = { head };

		while (true)
		{
			//list of reflected weights
			std::vector<weight> reflected_list;

			//iterate through master list
			for (auto i : master_list)
			{
				//iterate through selected simple roots that stabilize the weight
				for (size_t j = 0; j < simple.size(); j++)
				{
					//check if weight has ortho already, if not give it one
					if (i.ortho.size() == 0)
						i.ortho = to_ortho(i);

					/*
					* Generate reflected weight by reflecting i from masterlist by simple root j
					* Reflections are done in the ortho base, then we generate the omega base for
					* reflected weight (we neglect alpha for the moment)
					*/
					weight reflected;

					MatrixXr Reflector = ReflectionMatrix<rational<int>>(simple[j].ortho);
					reflected.ortho = Reflector * i.ortho;
					reflected.omega = to_omega(reflected);

					//Check to see if new reflected weight is unique
					auto findmaster = std::find(master_list.begin(), master_list.end(), reflected);
					auto findreflected = std::find(reflected_list.begin(), reflected_list.end(), reflected);
					if (findmaster == master_list.end() && findreflected == reflected_list.end())
						reflected_list.push_back(reflected);

				}
			}

			//If no new reflections are made, stop. All have been found
			if (reflected_list.empty())
				break;

			//concatenate reflected list with master list
			std::vector<weight> temp;
			temp.reserve(reflected_list.size() + master_list.size());
			temp.insert(temp.end(), master_list.begin(), master_list.end());
			temp.insert(temp.end(), reflected_list.begin(), reflected_list.end());
			master_list = temp;
		}
		return master_list;
	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::weylOrbit(weight head,std::vector<int> stabilizer)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "weylOrbit with stab called";
		std::vector<weight> master_list = { head };

		while (true)
		{
			//list of reflected weights
			std::vector<weight> reflected_list;

			//iterate through master list
			for (auto i : master_list)
			{
				//iterate through simple roots
				for (auto j : stabilizer)
				{
					//check if weight has ortho already, if not give it one
					if (i.ortho.size() == 0)
						i.ortho = to_ortho(i);

					/*
					* Generate reflected weight by reflecting i from masterlist by simple root j
					* Reflections are done in the ortho base, then we generate the omega base for
					* reflected weight (we neglect alpha for the moment)
					*/
					weight reflected;

					MatrixXr Reflector = ReflectionMatrix<rational<int>>(simple[j].ortho);
					reflected.ortho = Reflector * i.ortho;

					//Check to see if new reflected weight is unique
					auto findmaster = std::find(master_list.begin(), master_list.end(), reflected);
					auto findreflected = std::find(reflected_list.begin(), reflected_list.end(), reflected);
					if (findmaster == master_list.end() && findreflected == reflected_list.end())
						reflected_list.push_back(reflected);

				}
			}

			//If no new reflections are made, stop. All have been found
			if (reflected_list.empty())
				break;

			//concatenate reflected list with master list
			std::vector<weight> temp;
			temp.reserve(reflected_list.size() + master_list.size());
			temp.insert(temp.end(), master_list.begin(), master_list.end());
			temp.insert(temp.end(), reflected_list.begin(), reflected_list.end());
			master_list = temp;
		}
		return master_list;
	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::dominantWeights(std::vector<weight> representation)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "dominantWeights called";
		//Picks dominant weights from a representation and fills vector
		std::vector<weight> DominantWeights;

		for (auto i : representation)
		{
			//ensure omega basis
			if (i.omega.size() == 0)
				i.omega = to_omega(i);

			//dominant weight is omega basis with indicies all positive
			if (is_pos<rational<int>>(i.omega))
			{
				DominantWeights.push_back(i);

			}
		}
		return DominantWeights;
	}

	template<GroupType T>
	inline Matrix<rational<int>,Dynamic,1> LieBase<T>::to_alpha(weight V)
	{
		FILE_LOG(Diagnositics::logDEBUG4) << "to_alpha called";

		if (V.alpha.size() != 0)
			return V.alpha;
		else if (V.ortho.size() != 0)
			return ((PseudoInverse<rational<int>>(Cartan)).transpose())*((PseudoInverse<rational<int>>(Omega).transpose())*V.ortho);
		else
			return (PseudoInverse<rational<int>>(Cartan)).transpose() * V.omega;
	}

	template<GroupType T>
	inline Matrix<rational<int>, Dynamic, 1> LieBase<T>::to_omega(weight V)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "to_omega called";

		if (V.omega.size() != 0 )
			return V.omega;
		else if (V.ortho.size() != 0)
			return (PseudoInverse<rational<int>>(Omega).transpose())*V.ortho;
		else
			return Cartan.transpose()*V.alpha;
	}

	template<GroupType T>
	inline Matrix<rational<int>, Dynamic, 1> LieBase<T>::to_ortho(weight V)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "to_ortho called";

		if (V.ortho.size() != 0)
			return V.ortho;
		else if (V.omega.size() != 0)
			return Omega.transpose() * V.omega;
		else
			return Omega.transpose() * Cartan.transpose() * V.alpha;

	}

/*
*		Constructor
*/
	template<GroupType T>
	inline LieBase<T>::LieBase(const size_t Rank)
	{
		FILE_LOG(Diagnositics::logDEBUG2) << "LieBase Constructor rank "<<Rank<<" called";

		this->Rank = Rank;
		this->Group = T;
		createOrtho();
		createMatrices();
		createAllBases();

	}

	template<>
	inline void LieBase<GroupType::A>::createOrtho()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createOrtho A called";

		//Create Identity matrix
		MatrixXr Id_1 = Identity(Rank + 1);
		MatrixXr Id_2 = Identity(Rank);

		//Creates orthogonal simple roots, omega fundamental weights

		for (size_t i = 0; i < Rank; i++)
		{
			weight temp_r, temp_w;
			temp_r.ortho = Id_1.row(i) - Id_1.row(i + 1);
			temp_w.omega = Id_2.row(i);
			simple.push_back(temp_r);
			fweight.push_back(temp_w);
			FILE_LOG(Diagnositics::logDEBUG4) << temp_w.omega.transpose();

		}

		//Creates orthogonal positive roots
		//ei-ej
		for (size_t j = 1; j <= Rank; j++)
		{
			for (size_t i = 0; i < j; i++)
			{
				weight temp;
				temp.ortho = Id_1.row(i) - Id_1.row(j);
				positiver.push_back(temp);
			}
		}
	}

	template<>
	inline void LieBase<GroupType::B>::createOrtho()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createOrtho B called";

		//Create Identity matrix
		MatrixXr Id_1 = Identity(Rank);

		//Creates orthogonal simple roots, omega fundamental weights

		for (size_t i = 0; i < Rank; i++)
		{
			if (i < Rank - 1)
			{
				weight temp_r;
				temp_r.ortho = Id_1.row(i) - Id_1.row(i + 1);
				simple.push_back(temp_r);
			}
			else
			{
				weight temp_r;
				temp_r.ortho = Id_1.row(Rank - 1);
				simple.push_back(temp_r);
			}
			weight temp_w;
			temp_w.omega = Id_1.row(i);
			fweight.push_back(temp_w);
		}



		//Creates orthogonal positive roots
		// ei+ej,ei-ej,ei

		for (size_t i = 0; i < Rank; i++) {
			weight temp_p;
			temp_p.ortho = Id_1.row(i);
			positiver.push_back(temp_p);
		}

		for (size_t j = 1; j < Rank; j++) {
			for (size_t i = 0; i < j; i++) {
				weight temp1, temp2;
				temp1.ortho = Id_1.row(i) - Id_1.row(j);
				temp2.ortho = Id_1.row(i) + Id_1.row(j);
				positiver.push_back(temp1);
				positiver.push_back(temp2);
			}
		}

	}

	template<>
	inline void LieBase<GroupType::C>::createOrtho()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createOrtho C called";

		//Create Identity matrix
		MatrixXr Id_1 = Identity(Rank);

		//Creates orthogonal simple roots, omega fundamental weights

		for (size_t i = 0; i < Rank; i++)
		{
			if (i < Rank - 1)
			{
				weight temp_r;
				temp_r.ortho = Id_1.row(i) - Id_1.row(i + 1);
				simple.push_back(temp_r);
			}
			else
			{
				weight temp_r;
				temp_r.ortho = 2 * Id_1.row(Rank - 1);
				simple.push_back(temp_r);
			}
			weight temp_w;
			temp_w.omega = Id_1.row(i);
			fweight.push_back(temp_w);
		}

		//Creates orthogonal positive roots

		for (size_t i = 0; i < Rank; i++) {
			weight temp_p;
			temp_p.ortho = 2*Id_1.row(i);
			positiver.push_back(temp_p);
		}
		for (size_t j = 1; j < Rank; j++) {
			for (size_t i = 0; i < j; i++) {
				weight temp1, temp2;
				temp1.ortho = Id_1.row(i) - Id_1.row(j);
				temp2.ortho = Id_1.row(i) + Id_1.row(j);
				positiver.push_back(temp1);
				positiver.push_back(temp2);
			}
		}
	}

	template<>
	inline void LieBase<GroupType::D>::createOrtho()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createOrtho D called";

		//Create Identity matrix
		MatrixXr Id_1 = Identity(Rank);

		//Creates orthogonal simple roots, omega fundamental weights

		for (size_t i = 0; i < Rank; i++)
		{
			if (i < Rank - 1)
			{
				weight temp_r;
				temp_r.ortho = Id_1.row(i) - Id_1.row(i + 1);
				simple.push_back(temp_r);
			}
			else
			{
				weight temp_r;
				temp_r.ortho = Id_1.row(Rank - 1) + Id_1.row(Rank - 2);
				simple.push_back(temp_r);
			}
			weight temp_w;
			temp_w.omega = Id_1.row(i);
			fweight.push_back(temp_w);
		}

		//Creates orthogonal positive roots

		for (size_t i = 0; i < Rank; i++) {
			weight temp_p;
			temp_p.ortho = 2 * Id_1.row(i);
			positiver.push_back(temp_p);
		}
		for (size_t j = 1; j < Rank; j++) {
			for (size_t i = 0; i < j; i++) {
				weight temp1, temp2;
				if (j < Rank - 1)
				{
					temp2.ortho = Id_1.row(i) + Id_1.row(j);
					positiver.push_back(temp2);
				}

				temp1.ortho = Id_1.row(i) - Id_1.row(j);
				positiver.push_back(temp1);
			}
		}
	}

	template<>
	inline void LieBase<GroupType::G>::createOrtho()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createOrtho G called";

		//Creates orthogonal simple roots, 
		Matrix<rational<int>, 3, 1> simportho1(rational<int>(0),rational<int>(1),rational<int>(-1));
		Matrix<rational<int>, 3, 1> simportho2(rational<int>(1),rational<int>(-2),rational<int>(1));

		weight w1,w2;
		w1.ortho = simportho1;
		w2.ortho = simportho2;

		simple.push_back(w1);
		simple.push_back(w2);

		//Create Identity matrix
		MatrixXr Id_1 = Identity(Rank);

		//Creates omega fundamental weights
		for (size_t i = 0; i < Rank; i++)
		{
			weight temp_w;
			temp_w.omega = Id_1.row(i);
			fweight.push_back(temp_w);
		}

		//Creates orthogonal positive roots
		positiver = {w1,w2,w1+w2,2 * w1+w2,w1*3+w2};
		
	}

	template<>
	inline void LieBase<GroupType::F>::createOrtho()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createOrtho F called";

		//Creates orthogonal simple roots, 
		Matrix<rational<int>, 4, 1> simportho1(rational<int>(0),rational<int>(1),rational<int>(-1),rational<int>(0));
		Matrix<rational<int>, 4, 1> simportho2(rational<int>(0),rational<int>(0),rational<int>(1),rational<int>(-1));
		Matrix<rational<int>, 4, 1> simportho3(rational<int>(0),rational<int>(0),rational<int>(0),rational<int>(1));
		Matrix<rational<int>, 4, 1> simportho4(rational<int>(1,2),rational<int>(-1,2),rational<int>(-1,2),rational<int>(-1,2));

		weight w1,w2,w3,w4;
		w1.ortho = simportho1;
		w2.ortho = simportho2;
		w3.ortho = simportho3;
		w4.ortho = simportho4;

		simple.push_back(w1);
		simple.push_back(w2);
		simple.push_back(w3);
		simple.push_back(w4);

		//Create Identity matrix
		MatrixXr Id_1 = Identity(Rank);

		//Creates omega fundamental weights
		for (size_t i = 0; i < Rank; i++)
		{
			weight temp_w;
			temp_w.omega = Id_1.row(i);
			fweight.push_back(temp_w);
		}

		//Creates orthogonal positive roots

		//not implemented in this group
	}

	template<GroupType T>
	inline void LieBase<T>::createMatrices()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createMatrices called";

		Cartan.resize(Rank, Rank);
		//Cartan Defined
		for (size_t i = 0; i < Rank; i++)
		{
			for (size_t j = 0; j < Rank; j++)
			{
				Cartan(i, j) = master_formula<rational<int>>(simple[j].ortho, simple[i].ortho);
			}
		}

		CoCartan.resize(Rank, Rank);

		for (size_t i = 0; i < Rank; i++)
		{
			CoCartan.row(i) = 2 * simple[i].ortho / (simple[i].ortho.dot(simple[i].ortho));
		}
		Omega = (PseudoInverse<rational<int>>(CoCartan)).transpose();

		//QuadraticForm Defined
		QuadraticForm.resize(Rank, Rank);
		MatrixXr temp;
		temp.resize(Rank, Rank);
		for (size_t i = 0; i < Rank; i++)
		{
			temp(i, i) = rational<int>(1, 2)*simple[i].ortho.dot(simple[i].ortho);
		}

		QuadraticForm = RationalInverse<rational<int>>(Cartan) * temp;
	}

	template<>
	inline void LieBase<GroupType::A>::createMatrices()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createMatrices A called";

		Cartan.resize(Rank, Rank);
		//Cartan Defined
		for (size_t i = 0; i < Rank; i++)
		{
			for (size_t j = 0; j < Rank; j++)
			{
				Cartan(i, j) = master_formula<rational<int>>(simple[j].ortho, simple[i].ortho);
			}
		}

		CoCartan.resize(Rank, Rank + 1);

		for (size_t i = 0; i < Rank; i++)
		{
			CoCartan.row(i) = 2 * simple[i].ortho / (simple[i].ortho.dot(simple[i].ortho));
		}
		Omega = (PseudoInverse<rational<int>>(CoCartan)).transpose();

		//QuadraticForm Defined
		QuadraticForm.resize(Rank, Rank);
		MatrixXr temp;
		temp.resize(Rank, Rank);
		for (size_t i = 0; i < Rank; i++)
		{
			temp(i, i) = rational<int>(1, 2)*simple[i].ortho.dot(simple[i].ortho);
		}

		QuadraticForm = RationalInverse<rational<int>>(Cartan) * temp;
	}

	template<>
	inline void LieBase<GroupType::G>::createMatrices()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createMatrices G called";

		Cartan.resize(Rank, Rank);
		//Cartan Defined
		for (size_t i = 0; i < Rank; i++)
		{
			for (size_t j = 0; j < Rank; j++)
			{
				Cartan(i, j) = master_formula<rational<int>>(simple[j].ortho, simple[i].ortho);
			}
		}

		CoCartan.resize(Rank, Rank + 1);

		for (size_t i = 0; i < Rank; i++)
		{
			CoCartan.row(i) = 2 * simple[i].ortho / (simple[i].ortho.dot(simple[i].ortho));
		}
		Omega = (PseudoInverse<rational<int>>(CoCartan)).transpose();

		//QuadraticForm Defined
		QuadraticForm.resize(Rank, Rank);
		MatrixXr temp;
		temp.resize(Rank, Rank);
		for (size_t i = 0; i < Rank; i++)
		{
			temp(i, i) = rational<int>(1, 2)*simple[i].ortho.dot(simple[i].ortho);
		}

		QuadraticForm = RationalInverse<rational<int>>(Cartan) * temp;
	}

	template<GroupType T>
	inline void LieBase<T>::createAllBases()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "createAllBases called";

		//create all bases for simple roots
		for (size_t i = 0; i < simple.size(); i++)
		{
			simple[i].alpha = to_alpha(simple[i]);
			simple[i].omega = to_omega(simple[i]);
		}

		//create all bases for positive roots
		for (size_t i = 0; i < positiver.size(); i++)
		{
			positiver[i].alpha = to_alpha(positiver[i]);
			positiver[i].omega = to_omega(positiver[i]);
		}

		//Create rho, needed for dim and freudenthalRecursion
		//init rho vector = (1,1,1...)
		rho.resize(Rank);
		for (size_t i = 0; i < Rank; i++)
			rho(i) = rational<int>(1);

	}

	template<GroupType T>
	inline weight LieBase<T>::chamberRotate(weight w, int &counter)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "chamberRotate called";

		//ensure omega base and orthogonal base
		if (w.omega.size() == 0)
			w.omega = to_omega(w);
		if (w.ortho.size() == 0)
			w.ortho = to_ortho(w);

		//check to see if w is already in the dominant chamber
		if (is_pos<rational<int>>(w.omega))
			return w;

		weight reflected = w;

		//reflect to dominant chamber, counter flips +/- 1 each reflection
		while (true)
		{
			for (auto i : simple)
			{
				counter *= -1;

				//Reflection is done via ortho basis, then rotated to omega to check is_pos()
				weight temp;
				temp.ortho = ReflectionMatrix<rational<int>>(i.ortho) * reflected.ortho;
				temp.omega = to_omega(temp);
				reflected = temp;

				if (is_pos<rational<int>>(reflected.omega))
					return reflected;
			}
		}
		return w;
	}

	template<GroupType T>
	inline int LieBase<T>::freudenthalsRecursion(weight current, std::vector<weight> domWeights, std::vector<std::pair<int, weight>> stabilizedOrbits)
	{
		//Main implementation of the modified algorithm, only the stabilized orbits are included in summing
		FILE_LOG(Diagnositics::logDEBUG3) << "freudenthalsRecursion called";

		weight highest = domWeights[0];

		int k_current = k_lvl(current);
		int k_highest = k_lvl(highest);

		//check
		if (current == highest)
			return 1;

		int dummy = 0;
		int k = 1;
		rational<int> sum(1);

		//Recursion over the numerator
		do
		{
			for (auto stabs : stabilizedOrbits)
			{
				weight next;
				next.omega = current.omega + k * stabs.second.omega;
				weight rotated = chamberRotate(next,dummy);

				//Use chamberRotate() to rotated next to the dominant chamber, the counter isn't needed at this time
				auto in_rep = std::find(domWeights.begin(), domWeights.end(), rotated);

				if (in_rep != domWeights.end())
				{
					rational<int> inner_product = (rotated.omega.transpose() * QuadraticForm * stabs.second.omega);
					sum += inner_product * freudenthalsRecursion(rotated, domWeights, multiplicity(rotated)) * stabs.first;
					// TODO: Create static list of mulitplicites rather than recalculate previously visited weights
				}
			}
			k++;
		} while (k < k_highest - k_current);

		//Denominator Calculation
		rational<int> denom1 = (highest.omega + rho).transpose() * QuadraticForm * (highest.omega + rho);
		rational<int> denom2 = (current.omega + rho).transpose() * QuadraticForm * (current.omega + rho);

		rational<int> result = sum / (denom1 - denom2);

		int result_int = rational_cast<int>(result);

		return result_int;
	}

	template<GroupType T>
	inline std::vector<std::pair<int, weight>> LieBase<T>::multiplicity(weight highest)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "multiplicity called";

		/*
		* This is a helper function for the modified Freudenthal Recursion Formula
		* R. V. Moody, J. Patera, Fast Recursion Formula for Weight Multiplicities, Bull.Amer.Math.Soc.(N.S.)
			7 (1) (1982) 237�242.
		*/
		std::vector<std::pair<int, weight>> multiplicity;
		std::vector<int> stabilizers;

		//ensure omega
		if (highest.omega.size() == 0)
			highest.omega = to_omega(highest);

		//Stabilizer is where omega index is zero
		for (int i = 0; i < highest.omega.size(); i++)
			if (highest.omega(i) == 0)
				stabilizers.push_back(i);

		//Stabilized weights
		std::vector<weight> xis;

		//Loads xis (stabilized weights) with positive roots
		for (auto i : positiver)
		{
			bool xi_flag = true;

			for (auto j : stabilizers)
			{
				if (i.omega(j) >= 0 && xi_flag)
					xi_flag = true;
				else
					xi_flag = false;
			}
			if (xi_flag)
				xis.push_back(i);
		}

		//This consideres each xi, checks to see if its stabilizer is in root system based on simple roots labeled by 'stabilizer'
		//if so the overall multiplicity is equal to the size of the stabilized weylOrbit otherwise its a factor of 2
		for (auto i : xis)
		{
			std::vector<int> xi_stab;
			//This calculates the orbit size factor 1 or 2
			VectorXr temp = to_alpha(i);
			for (size_t d = 0; d < Rank; d++)
			{
				if (temp(d) > 0)
					xi_stab.push_back(d);
			}

			std::vector<int> set_dif(Rank);
			auto it = std::set_difference(xi_stab.begin(), xi_stab.end(), stabilizers.begin(), stabilizers.end(), set_dif.begin());
			set_dif.resize(it - set_dif.begin());

			if (set_dif.empty())
				multiplicity.push_back({ 1 * weylOrbit(i,stabilizers).size(), i });
			else
				multiplicity.push_back({ 2 * weylOrbit(i,stabilizers).size(), i });

		}


		return multiplicity;
	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::weightTower()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "weightTower() called";

		std::vector<weight> final_result;

		//concenate all simple orbits
		for (auto i : simple)
		{
			std::vector<weight> orbit = weylOrbit(i);
			final_result.insert(final_result.end(), orbit.begin(), orbit.end());
		}

		//this just gives a base sorting relation
		auto sort_comparator = [=](weight &X, weight &Y)->bool {
			if (X.omega.size() == 0)
				X.omega = to_omega(X);
			if (Y.omega.size() == 0)
				Y.omega = to_omega(Y);
			int i = 0;
			while (i < X.omega.size()) {
				if (X.omega(i) < Y.omega(i))
					return true;
				else if (X.omega(i) == Y.omega(i))
					i++;
				else
					return false;
			}
			return false;
		};

		std::sort(final_result.begin(), final_result.end(), sort_comparator);
		auto last = std::unique(final_result.begin(), final_result.end());
		final_result.erase(last, final_result.end());

		VectorXr z;
		z.resize(Rank);
		for (size_t i = 0; i < Rank; i++) z(i) = rational<int>(0);

		for (size_t i = 0; i < Rank; i++)
		{
			weight zeros;
			zeros.omega = z;
			final_result.push_back(zeros);
		}

		return final_result;

	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::weightTower(weight w)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "weightTower with weight called";

		//Given tower of weights for highest weight 'w'

		//Finding k levels -> think angular quantum numbers

		//init tower
		std::vector<weight> tower(2 *  k_lvl(w) + 1);
		tower[0] =  w ;

		//Build tower by subtracting simple roots from highest weight. Remove artifical degenerencies
		while (true)
		{
			std::vector<weight> temp;
			for (size_t i = 0; i < tower.size(); i++)
			{
				for (size_t j = 0; j < Rank; j++)
				{
					rational<int> p_val(1);
					while (p_val <= tower[i].omega(j))
					{
						//Check for duplicates
						weight current;
						current.omega = tower[i].omega - p_val * simple[j].omega;
						auto dup_chk1 = std::find(tower.begin(), tower.end(), current);
						auto dup_chk2 = std::find(temp.begin(), temp.end(), current);

						if (dup_chk1 == tower.end() && dup_chk2 == temp.end())
							tower.push_back(current);

						p_val++;
					}
				}
			}

			//exit condition -> we've reached the bottom
			if (temp.empty())
				break;
		}

		//We now need to use freudenthalRecursion to add proper degenerencies for dominantWeights in tower
		std::vector<weight> dom_weights = dominantWeights(tower);
		std::vector<std::pair<int, weight>> dom_weights_multi;


		//load list of degenerate weights
		for (auto wt : dom_weights)
		{
			int multi = freudenthalsRecursion(wt, dom_weights, multiplicity(wt));
			dom_weights_multi.push_back({ multi, wt });

		}

		//Modified algorithm states that tower is built from weylOrbits of dominantWeights with proper degenerencies
		std::vector<weight> final_result;

		for (auto i : weylOrbit(dom_weights_multi[0].second))
		{
			if (i.omega.size() == 0)
				i.omega = to_omega(i);
			std::cout << int_cast_Vector( i.omega).transpose() << std::endl;
		}
		for (auto i : dom_weights_multi)
		{
			for (int j = 0; j < i.first; j++)
			{
				std::vector<weight> orbit = weylOrbit(i.second);
				final_result.insert(final_result.end(), orbit.begin(), orbit.end());
			}
		}

		return final_result;
	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::tensorProductDecomp(weight w1, weight w2)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "tensorProductDecomp called";

		//Algorithm based on Klimyk's formula

		std::vector<weight> tower1 = weightTower(w1);
		std::vector<weight> mu;

		if (w1.omega.size() == 0)
			w1.omega = to_omega(w1);
		if (w2.omega.size() == 0)
			w2.omega = to_omega(w2);

		for (auto wt : tower1)
		{
			if (wt.omega.size() == 0)
				wt.omega = to_omega(wt);
			mu.push_back(weight(wt.omega + w2.omega + rho));
		}

		std::vector<std::pair<int, weight>> weight_parities;


		for (auto wt : mu)
		{
			int parity = 1;
			weight temp = chamberRotate(wt, parity);

			//ensure omega
			if (temp.omega.size() == 0)
				temp.omega = to_omega(temp);

			//Drop weights lying on chamber wall
			if (!has_zero(temp.omega))
			{
				temp.omega = temp.omega - rho;
				weight_parities.push_back({ parity, temp });
			}
		}

		//weight_partities are unsorted
		auto wp_sort = [](std::pair<int, weight> &X, std::pair<int, weight> &Y)->bool
		{
			int i = 0;
			while (i < X.second.omega.size())
			{
				if (X.second.omega(i) < Y.second.omega(i))
					return true;
				else if (X.second.omega(i) == Y.second.omega(i))
					i++;
				else
					return false;
			}
			return false;
		};
		std::sort(weight_parities.begin(), weight_parities.end(), wp_sort);

		auto multi_count = [](const std::pair<int, weight> &X, const std::pair<int, weight> &Y) { return X.second.omega == Y.second.omega; };
		auto it = std::adjacent_find(weight_parities.begin(), weight_parities.end(), multi_count);

		while (it != weight_parities.end())
		{
			int mult = it->first + (it + 1)->first;
			(it + 1)->first = mult;
			it->first = 0;
			it = std::adjacent_find(++it, weight_parities.end(), multi_count);
		}

		std::vector<weight> tensor_sum;
		for (auto i : weight_parities) {
			for (int j = 0; j < i.first; j++)
				tensor_sum.push_back(i.second);
		}

		for (auto i : tensor_sum)
			std::cout << " weight: " << i.omega.transpose() << std::endl;

		return tensor_sum;
	}

	template<GroupType T>
	size_t LieBase<T>::get_Rank()
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_Rank called";
		return Rank;
	}

	template<GroupType T>
	inline Matrix<rational<int>, Dynamic, Dynamic> LieBase<T>::get_Cartan() const
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_Cartan called";
		return Cartan;
	}

	template<GroupType T>
	inline Matrix<rational<int>, Dynamic, Dynamic> LieBase<T>::get_Omega() const
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_Omega called";
		return Omega;
	}

	template<GroupType T>
	inline Matrix<rational<int>, Dynamic, Dynamic> LieBase<T>::get_CoCartan() const
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_CoCartan called";
		return CoCartan;
	}

	template<GroupType T>
	inline Matrix<rational<int>, Dynamic, Dynamic> LieBase<T>::get_QuadraticForm() const
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_QuadraticForm called";
		return QuadraticForm;
	}

	template<GroupType T>
	inline int LieBase<T>::dim(weight w)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "dim called";

		//ensure omega basis
		if (w.omega.size() == 0)
			w.omega = to_omega(w);

		//main algorithm
		auto prod = rational<int>(1);
		for (auto i : positiver)
		{
			rational<int> numer = i.omega.transpose() * QuadraticForm * (w.omega + rho);
			rational<int> denom = rho.transpose() * QuadraticForm * (i.omega);
			prod *= numer / denom;
		}

		int final_ans = rational_cast<int>(prod);
		return final_ans;
	}

	template<GroupType T>
	inline int LieBase<T>::k_lvl(weight w)
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "k_lvl called";

		if (w.alpha.size() == 0)
			w.alpha = to_alpha(w);

		rational<int> k = w.alpha.sum();

		int final_ans = rational_cast<int>(k);
		return final_ans;
	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::get_simple() const
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_simple called";
		return this->simple;
	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::get_positiver() const
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_positiver called";
		return this->positiver;
	}

	template<GroupType T>
	inline std::vector<weight> LieBase<T>::get_fweight() const
	{
		FILE_LOG(Diagnositics::logDEBUG3) << "get_fweight called";
		return this->fweight;
	}

	template<>
	inline LieBase<GroupType::G>::LieBase(const size_t Rank)
	{
		FILE_LOG(Diagnositics::logDEBUG2) << "LieBase Constructor rank "<<Rank<<" called";
		
		if (Rank != 2)
		{
			FILE_LOG(Diagnositics::logWARNING) << "Lie Group 'G' is only defined for Rank 2. Rank was autoset 'Rank = 2'.";
		}

		this->Rank = 2;
		this->Group = GroupType::G;
		createOrtho();
		createMatrices();
		createAllBases();
	}

	template<>
	inline LieBase<GroupType::F>::LieBase(const size_t Rank)
	{
		FILE_LOG(Diagnositics::logDEBUG2) << "LieBase Constructor rank "<<Rank<<" called";
		
		if (Rank != 4)
		{
			FILE_LOG(Diagnositics::logWARNING) << "Lie Group 'F' is only defined for Rank 4. Rank was autoset 'Rank = 4'.";
		}

		this->Rank = 4;
		this->Group = GroupType::F;
		createOrtho();
		createMatrices();
		createAllBases();
	}

}



#endif // !REPRESENTATION_LIEALGEBRA_H
