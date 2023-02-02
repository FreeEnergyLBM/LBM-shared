#include <string>
#include <memory>
#include <utility>

class LBModel {
	 
public:
    template <typename T>
    LBModel(T&& obj): object(std::make_shared<Model<T>>(obj)){}
      
    double collide() const {
        return object->collide(); 
    }
    void precompute() const{
        object->precompute();
    }

    void initialise() const{
        object->initialise();
    }
	
   struct Concept {
       virtual ~Concept() {}
	   virtual double collide() const = 0;
       virtual void precompute() const = 0;
       virtual void initialise() const = 0;
   };

   template< typename T >
   struct Model : Concept {
       Model(const T& t) : object(t) {}
	   double collide() const override {
		   return object.collide();
	   }
       void precompute() const override {
		   object.precompute();
	   }
       void initialise() const override {
		   object.initialise();
	   }
     private:
       T object;
   };
private:
   std::shared_ptr<const Concept> object;
};

template<class stencil>
class CollisionBase{
    public:

        template<int i,int ...d>
        double computeGamma(const double *velocity) const;

        template<int... i>
        double computeFirstMoment(const double *distribution,std::index_sequence<i...> sequence) const;

        template<int d,int... i>
        double computeSecondMoment(const double *distribution,std::index_sequence<i...> sequence) const;

        double collideSRT(const double& old,const double& equilibrium,const double& tau) const;

        template<int i,int ...xyz>
        double forceSRT(const double& force,const double* velocity,const double& itau) const;

    private:
        template<int i,int ...d>
        double computeVelocityFactor(const double *velocity) const;

        enum{x=0,y=1,z=2};
        
};

template<class stencil>
template<int i,int ...d>
double CollisionBase<stencil>::computeGamma(const double *velocity) const{

    return stencil::ma_Weights[i]*(1.0+computeVelocityFactor<i,d...>(velocity));

};

template<class stencil>
template<int i,int ...xyz>
double CollisionBase<stencil>::computeVelocityFactor(const double *velocity) const{

    constexpr double ci_dot_velocity=((stencil::template ma_Cxyz<xyz>()[i]*velocity[xyz])+...);

    constexpr double velocity_dot_velocity=((velocity[xyz]*velocity[xyz])+...);

    return (ci_dot_velocity)/stencil::m_Cs2+(ci_dot_velocity*ci_dot_velocity)/(2.0*stencil::m_Cs2*stencil::m_Cs2)-(velocity_dot_velocity)/(2.0*stencil::m_Cs2);

};

template<class stencil>
template<int... i>
double CollisionBase<stencil>::computeFirstMoment(const double *distribution,std::index_sequence<i...> sequence) const{

    return ((distribution[i])+...);

}

template<class stencil>
template<int d,int... i>
double CollisionBase<stencil>::computeSecondMoment(const double *distribution,std::index_sequence<i...> sequence) const{

    return ((distribution[i]*stencil::template ma_Cxyz<d>[i])+...);

}

template<class stencil>
double CollisionBase<stencil>::collideSRT(const double& old,const double& equilibrium,const double& itau) const{

    return old-itau*(old-equilibrium);

}

template<class stencil>
template<int i,int ...xyz>
double CollisionBase<stencil>::forceSRT(const double& force,const double* velocity,const double& itau) const{

    constexpr double ci_dot_velocity=((stencil::template ma_Cxyz<xyz>()[i]*velocity[xyz])+...);
    
    return (1-DT*itau/2.0)*stencil::ma_Weights[i]*(((stencil::template ma_Cxyz<xyz>()[i]-velocity[xyz])/stencil::m_Cs2+ci_dot_velocity*stencil::template ma_Cxyz<xyz>()[i]/(stencil::m_Cs2*stencil::m_Cs2*force[xyz]))+...);

}