#include "RootStructure.h"

namespace Representation {
	using namespace Eigen;
	using namespace boost;	
	
	std::ostream & operator<<( std::ostream & out, const weight & w)
	{
		int n = w.omega.size();
		out << "(";
		for (int i = 0; i < n; i++)
		{
			out<< rational_cast<int>(w.omega(i));
			if (i < n - 1)
			{
				out << ",";
			}
		}
		out << ")";
		return out;
	}

	MatrixXi int_cast_Matrix(Matrix<rational<int>, Dynamic, Dynamic> X)
	{
		MatrixXi v = MatrixXi::Zero(X.size(),X.size());
		for (int i = 0; i < X.rows(); i++)
		{
			for (int j = 0; j < X.cols(); j++)
			{
				v(i,j) = rational_cast<int>(X(i,j));
			}
		}
		return v;
	}

	Matrix<rational<int>, Dynamic, Dynamic> Identity(int size)
	{
		{
			Matrix<rational<int>, Dynamic, Dynamic> M;
			M.resize(size, size);
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					if (i == j) M(i, j) = rational<int>(1);
					else M(i, j) = rational<int>(0);
				}
			}
			return M;
		}
	}

	VectorXi int_cast_Vector(Matrix<rational<int>, Dynamic, 1> X)
	{
		VectorXi v = VectorXi::Zero(X.size());
		for (int i = 0; i < X.size(); i++)
			v(i) = rational_cast<int>(X(i));

		return v;
	}

	bool weight::is_empty()
	{
		if (this->alpha.size() == 0 && this->omega.size()==0 && this->ortho.size()==0)
			return true;
		else
			return false;
	}

	bool weight::operator==(const weight & w)
	{
		if (this->omega.size() != 0 && w.omega.size() != 0)
		{
			if (this->omega == w.omega)
				return true;
			else
				return false;
		}		
		else if (this->alpha.size() != 0 && w.alpha.size() != 0)
		{
			if (this->alpha == w.alpha)
				return true;
			else
				return false;
		}
		else if (this->ortho.size() != 0 && w.ortho.size() != 0)
		{
			if (this->ortho == w.ortho)
				return true;
			else
				return false;
		}
		else
			return false;
	}

	weight & weight::operator=(const weight & w)
	{
		this->ortho = w.ortho;
		this->omega = w.omega;
		this->alpha = w.alpha;

		return *this;
	}

	//omega init weight
	weight::weight(VectorXr omega)
	{
		this->omega = omega;
	}
	

}