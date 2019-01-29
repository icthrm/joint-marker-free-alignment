#include "opts.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "mrcmx/mrcstack.h"

#define MILLION 1000000

#define PI_180 0.01745329252f

#ifndef PI
#define     PI  3.14159265358979323846
#endif

#define D2R(__ANGLE__) ((__ANGLE__) * PI_180)


struct Coeff{
	union{
		double p[20];
		struct{
			double a[10];
			double b[10];
		};
	};
	
};

bool ReadAngles(std::vector<float>& angles, const char* name)
{
    std::ifstream in(name);
    if(!in.good()){
        return false;
    }

    while(in.good()){
		float val;
		in>>val;
		if(in.fail()) {
			break;
		}
		
        angles.push_back(val);
    }
    in.close();
    return true;
}

void TranslateAngleToCoefficients(const std::vector<float>& angles, const std::vector<float>& xangles, std::vector<Coeff>& coeffs){
	coeffs.resize(angles.size());
	for(int i = 0; i < angles.size(); i++){
		memset(coeffs[i].p, 0, sizeof(double)*20);
		float beta = D2R(angles[i]);
		float alpha = D2R(xangles[i]);
		
		coeffs[i].a[0] = 0; //
		coeffs[i].a[1] = cos(beta); //x
		coeffs[i].a[2] = sin(alpha)*sin(beta); //y
		coeffs[i].a[3] = -cos(alpha)*sin(beta); //z
		coeffs[i].b[0] = 0; //
		coeffs[i].b[1] = 0; //x
		coeffs[i].b[2] = cos(alpha); //y
		coeffs[i].b[3] = sin(alpha); //z
	}
}

/*solve inverse transfroms defined by Geometry; substitute the inversion into coefficients*/
void DecorateCoefficients(std::vector<Coeff>& coeffs, const Geometry& geo)
{
	double alpha = -D2R(geo.pitch_angle), beta = D2R(geo.offset), t = -geo.zshift;
	double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
	double ca2 = ca*ca, sa2 = sa*sa, cb2 = cb*cb, sb2 = sb*sb;
	
	for(int i = 0; i < coeffs.size(); i++){
		double a[10], b[10];
		memcpy(a, coeffs[i].a, sizeof(double)*10);
		memcpy(b, coeffs[i].b, sizeof(double)*10);
		coeffs[i].a[0] = a[0];
		coeffs[i].a[1] = (a[2]*sa*sb+a[3]*ca*sb+a[1]*cb);//*x
		coeffs[i].a[2] = (a[2]*ca-a[3]*sa);//*y
		coeffs[i].a[3] = (a[3]*ca*cb+a[2]*cb*sa-a[1]*sb);//*z
		coeffs[i].a[4] = (a[4]*ca*cb-a[5]*cb*sa-a[6]*sa2*sb-2*a[9]*ca*sa*sb+a[6]*ca2*sb+2*a[8]*ca*sa*sb);//*x*y
		coeffs[i].a[5] = (a[4]*cb2*sa-a[5]*ca*sb2+a[5]*ca*cb2-2*a[7]*cb*sb-a[4]*sa*sb2+2*a[9]*ca2*cb*sb+2*a[8]*cb*sa2*sb+2*a[6]*ca*cb*sa*sb);//*x*z
		coeffs[i].a[6] = (2*a[8]*ca*cb*sa-a[4]*ca*sb+a[5]*sa*sb+a[6]*ca2*cb-a[6]*cb*sa2-2*a[9]*ca*cb*sa);//*y*z
		coeffs[i].a[7] = (a[7]*cb2+a[5]*ca*cb*sb+a[9]*ca2*sb2+a[8]*sa2*sb2+a[4]*cb*sa*sb+a[6]*ca*sa*sb2);//*x^2
		coeffs[i].a[8] = (a[8]*ca2+a[9]*sa2-a[6]*ca*sa);//*y^2
		coeffs[i].a[9] = (a[7]*sb2+a[9]*ca2*cb2+a[8]*cb2*sa2-a[4]*cb*sa*sb+a[6]*ca*cb2*sa-a[5]*ca*cb*sb);//*z^2
		
		coeffs[i].b[0] = b[0];
		coeffs[i].b[1] = (b[2]*sa*sb+b[3]*ca*sb+b[1]*cb);//*x
		coeffs[i].b[2] = (b[2]*ca-b[3]*sa);//*y
		coeffs[i].b[3] = (b[3]*ca*cb+b[2]*cb*sa-b[1]*sb);//*z
		coeffs[i].b[4] = (b[4]*ca*cb-b[5]*cb*sa-b[6]*sa2*sb-2*b[9]*ca*sa*sb+b[6]*ca2*sb+2*b[8]*ca*sa*sb);//*x*y
		coeffs[i].b[5] = (b[4]*cb2*sa-b[5]*ca*sb2+b[5]*ca*cb2-2*b[7]*cb*sb-b[4]*sa*sb2+2*b[9]*ca2*cb*sb+2*b[8]*cb*sa2*sb+2*b[6]*ca*cb*sa*sb);//*x*z
		coeffs[i].b[6] = (2*b[8]*ca*cb*sa-b[4]*ca*sb+b[5]*sa*sb+b[6]*ca2*cb-b[6]*cb*sa2-2*b[9]*ca*cb*sa);//*y*z
		coeffs[i].b[7] = (b[7]*cb2+b[5]*ca*cb*sb+b[9]*ca2*sb2+b[8]*sa2*sb2+b[4]*cb*sa*sb+b[6]*ca*sa*sb2);//*x^2
		coeffs[i].b[8] = (b[8]*ca2+b[9]*sa2-b[6]*ca*sa);//*y^2
		coeffs[i].b[9] = (b[7]*sb2+b[9]*ca2*cb2+b[8]*cb2*sa2-b[4]*cb*sa*sb+b[6]*ca*cb2*sa-b[5]*ca*cb*sb);//*z^2
	}
	
	//considering z_shift
	for(int i = 0; i < coeffs.size(); i++){
		double a[10], b[10];
		memcpy(a, coeffs[i].a, sizeof(double)*10);
		memcpy(b, coeffs[i].b, sizeof(double)*10);
		
		coeffs[i].a[0] = a[0]+a[3]*t+a[9]*t*t;
		coeffs[i].a[1] = a[1]+a[5]*t;//*x
		coeffs[i].a[2] = a[2]+a[6]*t;//*y
		coeffs[i].a[3] = a[3]+2*a[9]*t;//*z
		
		coeffs[i].b[0] = b[0]+b[3]*t+b[9]*t*t;
		coeffs[i].b[1] = b[1]+b[5]*t;//*x
		coeffs[i].b[2] = b[2]+b[6]*t;//*y
		coeffs[i].b[3] = b[3]+2*b[9]*t;//*z
	}
}

void (*funL)(const Coeff&, double, double, double, double*);

void WarpPosition(const Coeff& coeff, double X, double Y, double Z, double* n){
	n[0] = coeff.a[0]+coeff.a[1]*X+coeff.a[2]*Y+coeff.a[3]*Z+coeff.a[4]*X*Y+coeff.a[5]*X*Z+coeff.a[6]*Y*Z+coeff.a[7]*X*X+coeff.a[8]*Y*Y+coeff.a[9]*Z*Z;
	n[1] = coeff.b[0]+coeff.b[1]*X+coeff.b[2]*Y+coeff.b[3]*Z+coeff.b[4]*X*Y+coeff.b[5]*X*Z+coeff.b[6]*Y*Z+coeff.b[7]*X*X+coeff.b[8]*Y*Y+coeff.b[9]*Z*Z;
}

void LinearPosition(const Coeff& coeff, double X, double Y, double Z, double* n){
	n[0] = coeff.a[0]+coeff.a[1]*X+coeff.a[2]*Y+coeff.a[3]*Z;
	n[1] = coeff.b[0]+coeff.b[1]*X+coeff.b[2]*Y+coeff.b[3]*Z;
}

void ValCoef(const Point3DF& origin, const Point3D& coord, const Coeff& coeff, Weight* wt)
{
	double x, y;

	double X, Y, Z, n[2];
	X = coord.x-origin.x; Y = coord.y-origin.y; Z = coord.z-origin.z;
	
// 	funL(coeff, X, Y, Z, n);
	n[0] = coeff.a[0]+coeff.a[1]*X+coeff.a[2]*Y+coeff.a[3]*Z;
	n[1] = coeff.b[0]+coeff.b[1]*X+coeff.b[2]*Y+coeff.b[3]*Z;
	
	x = n[0]+origin.x; y = n[1]+origin.y;

	wt->x_min = floor(x);
	wt->y_min = floor(y);

	wt->x_min_del = x - wt->x_min;
	wt->y_min_del = y - wt->y_min;
}

void Reproject(const Point3DF& origin, const Volume& vol, const Coeff& coeff, Slice& reproj_val, Slice& reproj_wt){

	Point3D coord; int n;

	for(int z = 0; z < vol.height; z++){
		float* vdrefz = vol.data+z*(size_t)vol.width*vol.length;
		coord.z = z+vol.z;
		
		for(int y = 0; y < vol.length; y++){
			float* vdrefy = vdrefz+y*vol.width;
			coord.y = y+vol.y;
			
			for(int x = 0; x < vol.width; x++){
				float* vdrefx = vdrefy+x;
				coord.x = x+vol.x; Weight wt;
				ValCoef(origin, coord, coeff, &wt);

				if(wt.x_min >= 0 && wt.x_min < vol.width && wt.y_min >= 0 && wt.y_min < vol.length){ //(x_min, y_min)
					n = wt.x_min + wt.y_min * vol.width; //index in reproj
					reproj_val.data[n] += (1-wt.x_min_del) * (1-wt.y_min_del) * (*vdrefx);
					reproj_wt.data[n] += (1-wt.x_min_del) * (1-wt.y_min_del);
				}
				if((wt.x_min+1) >= 0 && (wt.x_min+1) < vol.width && wt.y_min >= 0 && wt.y_min < vol.length){ //(x_min+1, y_min)
					n = wt.x_min+1 + wt.y_min * vol.width; //index in reproj
					reproj_val.data[n] += wt.x_min_del * (1-wt.y_min_del) * (*vdrefx);
					reproj_wt.data[n] += wt.x_min_del * (1-wt.y_min_del);
				}
				if(wt.x_min >= 0 && wt.x_min < vol.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < vol.length){ //(x_min, y_min+1)
					n = wt.x_min + (wt.y_min+1) * vol.width; //index in reproj
					reproj_val.data[n] += (1-wt.x_min_del) * wt.y_min_del * (*vdrefx);
					reproj_wt.data[n] += (1-wt.x_min_del) * wt.y_min_del;
				}
				if((wt.x_min+1) >= 0 && (wt.x_min+1) < vol.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < vol.length){ //(x_min+1, y_min+1)
					n = (wt.x_min+1) + (wt.y_min+1) * vol.width; //index in reproj
					reproj_val.data[n] += wt.x_min_del * wt.y_min_del * (*vdrefx);
					reproj_wt.data[n] += wt.x_min_del * wt.y_min_del;
				}
			}
		}
	}
}

inline void BilinearValue(const Slice& slc, const Weight& wt, float* val, float* vwt)
{
	int n;
	if(wt.x_min >= 0 && wt.x_min < slc.width && wt.y_min >= 0 && wt.y_min < slc.height){ //(x_min, y_min)
		n = wt.x_min + wt.y_min * slc.width;
		*val += (1-wt.x_min_del) * (1-wt.y_min_del) * slc.data[n];
		*vwt += (1-wt.x_min_del) * (1-wt.y_min_del);
	}
	if((wt.x_min+1) >= 0 && (wt.x_min+1) < slc.width && wt.y_min >= 0 && wt.y_min < slc.height){ //(x_min+1, y_min)
		n = wt.x_min+1 + wt.y_min * slc.width;
		*val += wt.x_min_del * (1-wt.y_min_del) * slc.data[n];
		*vwt += wt.x_min_del * (1-wt.y_min_del);
	}
	if(wt.x_min >= 0 && wt.x_min < slc.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < slc.height){ //(x_min, y_min+1)
		n = wt.x_min + (wt.y_min+1) * slc.width;
		*val += (1-wt.x_min_del) * wt.y_min_del * slc.data[n];
		*vwt += (1-wt.x_min_del) * wt.y_min_del;
	}
	if((wt.x_min+1) >= 0 && (wt.x_min+1) < slc.width && (wt.y_min+1) >= 0 && (wt.y_min+1) < slc.height){ //(x_min+1, y_min+1)
		n = wt.x_min+1 + (wt.y_min+1) * slc.width;
		*val += wt.x_min_del * wt.y_min_del * slc.data[n];
		*vwt += wt.x_min_del * wt.y_min_del;
	}
}

void BackProject(const Point3DF& origin, MrcStackM& projs, Volume& vol, Coeff coeffv[])
{
	Slice proj(projs.X(), projs.Y()); 
	Point3D coord;

	memset(vol.data, 0, sizeof(float)*vol.length*vol.width*vol.height);
	
	for(int idx = 0; idx < projs.Z(); idx++){
		printf("BPT begin to read %d projection for %d z-coordinate\n", idx, vol.z);
		projs.ReadSliceZ(idx, proj.data);

		for(int z = 0; z < vol.height; z++){
			float* vdrefz = vol.data+z*(size_t)vol.width*vol.length;
			coord.z = z+vol.z;
			
			for(int y = 0; y < vol.length; y++){
				float* vdrefy = vdrefz+y*vol.width;
				coord.y = y+vol.y;
				
				for(int x = 0; x < vol.width; x++){
					coord.x = x+vol.x; 
					Weight wt; float s = 0, c = 0;
					
					ValCoef(origin, coord, coeffv[idx], &wt);
					BilinearValue(proj, wt, &s, &c);          

					if(c){
						*(vdrefy+x) += (float)(s / c);
					}
				}
			}
		}
	}
} 

void UpdateVolumeByProjDiff(const Point3DF& origin, const Slice& diff, Volume& vol, float gamma, const Coeff& coeff){
	
	Point3D coord;

	for(int z = 0; z < vol.height; z++){
		float* vdrefz = vol.data+z*(size_t)vol.width*vol.length;
		coord.z = z+vol.z;
		
		for(int y = 0; y < vol.length; y++){
			float* vdrefy = vdrefz+y*vol.width;
			coord.y = y+vol.y;
			
			for(int x = 0; x < vol.width; x++){
				coord.x = x+vol.x;
				Weight wt; float s = 0, c = 0;
					
				ValCoef(origin, coord, coeff, &wt);
				BilinearValue(diff, wt, &s, &c);   

				if(c){
					*(vdrefy+x) += (float) (s / c)*gamma;
				}
			}
		}
	}
}

void UpdateWeightsByProjDiff(const Point3DF& origin, const Slice& diff, Volume& values, Volume& weights, const Coeff& coeff){
	
	Point3D coord;

	for(int z = 0; z < values.height; z++){
		float* vdrefz = values.data+z*(size_t)values.width*values.length;
		float* wdrefz = weights.data+z*(size_t)weights.width*weights.length;
		coord.z = z+values.z;
		
		for(int y = 0; y < values.length; y++){
			float* vdrefy = vdrefz+y*values.width;
			float* wdrefy = wdrefz+y*weights.width;
			coord.y = y+values.y;
			
			for(int x = 0; x < values.width; x++){
				coord.x = x+values.x;
				Weight wt; float s = 0, c = 0;
					
				ValCoef(origin, coord, coeff, &wt);
				BilinearValue(diff, wt, &s, &c);   

				*(vdrefy+x) += s;
				*(wdrefy+x) += c;
			}
		}
	}
}

void UpdateVolumeByWeights(Volume& vol, Volume& values, Volume& weights, float gamma){
	size_t pxsize = vol.height*vol.width*vol.length;
	
	for(size_t i = pxsize; i--;){
		if(weights.data[i]){
			vol.data[i] += values.data[i]/weights.data[i]*gamma;
		}
	}
}

void SART(const Point3DF& origin, MrcStackM& projs, Volume& vol, Coeff coeffv[], int iteration, float gamma)
{
	int pxsize = projs.X()*projs.Y();
	Slice reproj_val(projs.X(), projs.Y());	//reprojection value
	Slice reproj_wt(projs.X(), projs.Y());	//reprojection weight
	Slice projection(projs.X(), projs.Y());
	
	for(int i = 0; i < iteration; i++){
		for(int idx = 0; idx < projs.Z(); idx++){
			memset(reproj_val.data, 0, sizeof(float)*pxsize);
			memset(reproj_wt.data, 0, sizeof(float)*pxsize);
			
			Reproject(origin, vol, coeffv[idx], reproj_val, reproj_wt);
		
			MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			
			printf("SART begin to read %d projection for %d z-coordinate (iteraion %d)\n", idx, vol.z, i);
			projs.ReadSliceZ(idx, projection.data);
			
			for(int n = 0; n < pxsize; n++){
				if(reproj_wt.data[n] != 0){
					reproj_val.data[n] /= reproj_wt.data[n];
				}
				reproj_val.data[n] = projection.data[n]-reproj_val.data[n];
			}
			
			UpdateVolumeByProjDiff(origin, reproj_val, vol, gamma, coeffv[idx]);
		}
	}
}

void SIRT(const Point3DF& origin, MrcStackM& projs, Volume& vol, Coeff coeffv[], int iteration, float gamma)
{
	int pxsize = projs.X()*projs.Y();
	Slice reproj_val(projs.X(), projs.Y());	//reprojection value
	Slice reproj_wt(projs.X(), projs.Y());	//reprojection weight
	Slice projection(projs.X(), projs.Y());
	
	Volume valvol(vol.x, vol.y, vol.z, vol.length, vol.width, vol.height);
	Volume wtvol(vol.x, vol.y, vol.z, vol.length, vol.width, vol.height);
	
	for(int i = 0; i < iteration; i++){
		memset(valvol.data, 0, sizeof(float)*valvol.length*valvol.width*valvol.height);	
		memset(wtvol.data, 0, sizeof(float)*wtvol.length*wtvol.width*wtvol.height);	
		
		for(int idx = 0; idx < projs.Z(); idx++){
			memset(reproj_val.data, 0, sizeof(float)*pxsize);
			memset(reproj_wt.data, 0, sizeof(float)*pxsize);
			
			Reproject(origin, vol, coeffv[idx], reproj_val, reproj_wt);		//vol is not changed during iteration
		
			MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			
			printf("SIRT begin to read %d projection for %d z-coordinate (iteraion %d)\n", idx, vol.z, i);
			projs.ReadSliceZ(idx, projection.data);
			
			for(int n = 0; n < pxsize; n++){
				if(reproj_wt.data[n]){
					reproj_val.data[n] /= reproj_wt.data[n];
				}
				reproj_val.data[n] = projection.data[n]-reproj_val.data[n];
			}
			
			UpdateWeightsByProjDiff(origin, reproj_val, valvol, wtvol, coeffv[idx]);
		}
		
		UpdateVolumeByWeights(vol, valvol, wtvol, gamma);
	}
}

struct SysInfo{
	int id;
	int procs;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int namelen;
};

int ATOM(options& opt, int myid, int procs){
	
	MrcStackM projs, mrcvol;
	if(!projs.ReadFile(opt.input)){
		printf("File %s cannot access.\n", opt.input);
		
		return -1;
	}
	
	if(myid == 0){
		projs.ReadHeader();
	}
	MPI_Bcast(&(projs.header), sizeof(MRCheader), MPI_CHAR, 0, MPI_COMM_WORLD);
	
	mrcvol.InitializeHeader();
	mrcvol.SetSize(projs.X(), projs.Y(), opt.thickness);

	std::vector<float> angles;
	ReadAngles(angles, opt.angle);
	
	std::vector<float> xangles;
	
	if(opt.xangle[0] != '\0'){
		ReadAngles(xangles, opt.xangle);
	}
	else{
		xangles.resize(angles.size(), 0.0);
	}
	
	std::vector<Coeff> params;
	TranslateAngleToCoefficients(angles, xangles, params);
	
	Geometry geo;
	geo.offset = opt.offset;
	geo.pitch_angle = opt.pitch_angle;
	geo.zshift = opt.zshift;
	
	DecorateCoefficients(params, geo);
	
	int height;
	int zrem = mrcvol.Z()%procs;
	int volz;   //the start slice of reproject per process
	
	if(myid < zrem){
		height = mrcvol.Z()/procs+1;
		volz = height * myid;
	}
	else{
		height = mrcvol.Z()/procs;
		volz = height * myid+zrem;
	}

	Volume vol(0, 0, volz, mrcvol.Y(), mrcvol.X(), height);

	std::cout<<myid<<": ("<<vol.x<<","<<vol.y<<","<<vol.z<<")"<<"&("<<vol.width<<","<<vol.length<<","<<vol.height<<")"<<std::endl;

	Point3DF origin;
	
	origin.x = mrcvol.X()*.5;
	origin.y = mrcvol.Y()*.5;
	origin.z = mrcvol.Z()*.5;

	if(myid == 0){
		printf("origin.x is %f, origin.y is %f, origin.z is %f\n", origin.x, origin.y, origin.z);
	}
	
	/***************************reprojection************************/
	if(opt.method == "RP"){
		projs.ReadBlock(vol.z, vol.z+vol.height, 'z', vol.data);
		
		mrcvol.SetSize(projs.X(), projs.Y(), params.size());
		if(opt.output[0] == '\0'){
			const char* name = "reproj.mrc";
			memcpy(opt.output, name, sizeof(char)*11);
		}
		
		mrcvol.WriteToFile(opt.output);
	
		if(myid == 0){
			mrcvol.WriteHeader();
		}
		
		int pxsize = projs.X()*projs.Y();
		Slice reproj_val(projs.X(), projs.Y());	//reprojection value
		Slice reproj_wt(projs.X(), projs.Y());	//reprojection weight
			
		for(int idx = 0; idx < params.size(); idx++){
			memset(reproj_val.data, 0, sizeof(float)*pxsize);
			memset(reproj_wt.data, 0, sizeof(float)*pxsize);
			
			Reproject(origin, vol, params[idx], reproj_val, reproj_wt);
			MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
			
			for(int n = 0; n < pxsize; n++){
				if(reproj_wt.data[n] != 0){
					reproj_val.data[n] /= reproj_wt.data[n];
				}
			}
			
			if(myid == 0){
				printf("RP begin to write %d projection\n", idx);
				mrcvol.WriteSlice(idx, 'z', reproj_val.data);
			}
		}
		
		if(myid == 0){
			mrcvol.UpdateHeader();
		}

		projs.Close();
		mrcvol.Close();

		return 0;
	}

	/**********************reconstruction along Z-axis*******************/
	if(opt.method == "BPT"){
		if(opt.output[0] == '\0'){
			const char* name = "bpt.rec";
			memcpy(opt.output, name, sizeof(char)*8);
		}
		
		BackProject(origin, projs, vol, &params[0]);
	}
	else if(opt.method == "SART"){
		const char* reference = opt.initial;
		if(reference[0] == '\0'){
			memset(vol.data, 0, vol.width*vol.length*vol.height*sizeof(float));
		}
		else{
			MrcStackM init;
			init.ReadFile(reference);
			init.ReadHeader();
			init.ReadBlock(vol.z, vol.z+vol.height, 'z', vol.data);
			init.Close();
		}
		
		if(opt.output[0] == '\0'){
			const char* name = "sart.rec";
			memcpy(opt.output, name, sizeof(char)*9);
		}
		
		SART(origin, projs, vol, &params[0], opt.iteration, opt.gamma);
	}
	else if(opt.method == "SIRT"){
		const char* reference = opt.initial;
		if(reference[0] == '\0'){
			memset(vol.data, 0, vol.width*vol.length*vol.height*sizeof(float));
		}
		else{
			MrcStackM init;
			init.ReadFile(reference);
			init.ReadHeader();
			init.ReadBlock(vol.z, vol.z+vol.height, 'z', vol.data);
			init.Close();
		}

		if(opt.output[0] == '\0'){
			const char* name = "sirt.rec";
			memcpy(opt.output, name, sizeof(char)*9);
		}
		
		SIRT(origin, projs, vol, &params[0], opt.iteration, opt.gamma);
	}

	mrcvol.WriteToFile(opt.output);
	
	if(myid == 0){
		mrcvol.WriteHeader();
	}

	mrcvol.WriteBlock(vol.z, vol.z+height, 'z', vol.data);

	MPI_Barrier(MPI_COMM_WORLD);
	
	if(myid == 0){
		mrcvol.UpdateHeader();
	}

	projs.Close();
	mrcvol.Close();

	return 0;
}

int main(int argc, char *argv[]){
	SysInfo info;

	MPI_Init(&argc, &argv); //parallel init
	MPI_Comm_rank(MPI_COMM_WORLD, &(info.id));
	MPI_Comm_size(MPI_COMM_WORLD, &(info.procs));
	MPI_Get_processor_name(info.processor_name, &(info.namelen));

	options opts;
    InitOpts(&opts);

    if(GetOpts(argc, argv, &opts) <= 0) {
        EX_TRACE("***WRONG INPUT.\n");
        return -1;
    }
    
    if(info.id == 0){
		PrintOpts(opts);
	}

	ATOM(opts, info.id, info.procs);
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();                 //parallel finish

    return 0;
}

