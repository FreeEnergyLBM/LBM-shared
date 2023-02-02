
template<class stencil>
class BodyForce{
    public:
        double compute() const;
        void precompute() const;
        double computeDensitySource() const;
        double computeVelocitySource() const;
    private:
        double magnitude;
};

template<class stencil>
double BodyForce<stencil>::compute() const{
    return magnitude*Density;
}

template<class stencil>
void BodyForce<stencil>::precompute() const{

}

template<class stencil>
double BodyForce<stencil>::computeDensitySource() const{
    
}

template<class stencil>
double BodyForce<stencil>::computeVelocitySource() const{
    
}
