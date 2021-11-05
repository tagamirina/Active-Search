#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define FNL 2048
#define XY1D(x, y, w)	((x) + (y) * w)
#define ip_safe(x, y, w, h) 	((x) >= 0 && (x) < (w) && (y) >= 0 && (y) < (h))
#define rint(x)		((x >= 0)?((int)(x + 0.5)):((int)(x - 0.5)))

double calc_match_rate(cv::Mat *S_color, cv::Mat *II_color, int j, int i, int skip_n, unsigned char *isrc, unsigned char *itmp, int t_leng, int s_leng) ;
double calc_match_hist(cv::Mat *S_color, cv::Mat *II_color, int j, int i, int skip_n, unsigned char *isrc, unsigned char *itmp, int t_leng, int s_leng) ;
void NcorrBy1D( cv::Mat *I, cv::Mat *S, int si, int sj ) ;

int main (int argc,char **argv){

	int i, j, idx, inum, points, count, found, fmaxx ;
	cv::Mat II_color, SS_color, actmap, I_color, S_color ;
	cv::Mat	RESR, RESG, RESB ;
	char ifname[FNL], iifname[FNL], resfname[FNL], ssfname[FNL], mapfname[FNL] ;
	double uplimit, upl_mag, ref_r, tgt_r, rate, s_rate, image_mag ;
	int isx, isy ,ssx, ssy, skips_n, ref_s, n, nn ;
	int sch_x1, sch_y1, sch_x2, sch_y2 ;
	int tgt_x1, tgt_y1, tgt_x2, tgt_y2, tgt_s ;
	int com_x1, com_y1, com_x2, com_y2, com_s, dif_s ;
	int si, sj, ei, ej ;
	unsigned char *isrc, *itmp ;
	int t_leng, s_leng ;
	int ref_x1, ref_x2, ref_y1, ref_y2 ;
	FILE *fp ;

	skips_n = 1; // Width to skip
	count = 0 ;
	upl_mag = 0.75 ; // Multiplier for active search limit
	image_mag = 1.0 ; // Image size reduction rate
	ref_r = tgt_r = 0.0 ;
	si = sj = 0 ;
	s_rate = 0.0 ;
	found = 0 ;
	fmaxx = 1024 ;
	n = 0 ;
	nn = 0 ;

	for( idx = 1 ; idx < argc ; idx++ ){
		if(      !strcmp( argv[idx], "-i" )) strcpy( ifname, argv[++idx] ) ;
		if(      !strcmp( argv[idx], "-ii" )) strcpy( iifname, argv[++idx] ) ;
		else if( !strcmp( argv[idx], "-inum" )) inum = atoi( argv[++idx] ) ;
		else if( !strcmp( argv[idx], "-points" )) points = atoi( argv[++idx] ) ;
		else if( !strcmp( argv[idx], "-ss" )) strcpy( ssfname, argv[++idx] ) ;
		else if( !strcmp( argv[idx], "-res" )) strcpy( resfname, argv[++idx] ) ;
		else if( !strcmp( argv[idx], "-map" )) strcpy( mapfname, argv[++idx] ) ;
	}

	II_color = cv::imread( ifname, cv::IMREAD_COLOR ) ;
	SS_color = cv::imread( iifname, cv::IMREAD_COLOR ) ;
	isx = II_color.cols ;
	isy = II_color.rows ;
	ssx = SS_color.cols ;
	ssy = SS_color.rows ;
	sch_x1 = -isx / 2 ;
	sch_y1 = -isy / 2 ;
	sch_x2 = ssx - (isx / 2) ;
	sch_y2 = ssy - (isy / 2) ;
	actmap = (cv::Mat_<cv::Vec3b>(ssy, ssx)) ;
	RESR = (cv::Mat_<unsigned char>(ssy, ssx));
	RESG = (cv::Mat_<unsigned char>(ssy, ssx));
	RESB = (cv::Mat_<unsigned char>(ssy, ssx));


	// Calculating the width of a line
	t_leng = isx*3;
	s_leng = ssx*3;
	isrc = (unsigned char*)malloc(s_leng * ssy * sizeof(unsigned char)) ;
	itmp = (unsigned char*)malloc(t_leng * isy * sizeof(unsigned char)) ;

    // Initialization
	for(j = 0;j < ssy;j++){
		for(i = 0;i < ssx;i++){
			isrc[n] = SS_color.at<cv::Vec3b>(j, i)[0] ;
			isrc[n + 1] = SS_color.at<cv::Vec3b>(j, i)[1] ;
			isrc[n + 2] = SS_color.at<cv::Vec3b>(j, i)[2] ;
			n += 3 ;
		}
	}
	//printf("n %d\n",n) ;
	n = 0 ;
	for(j = 0;j < isy;j++){
		for(i = 0;i < isx;i++){
			itmp[n] = II_color.at<cv::Vec3b>(j, i)[0] ;
			itmp[n + 1] = II_color.at<cv::Vec3b>(j, i)[1] ;
			itmp[n + 2] = II_color.at<cv::Vec3b>(j, i)[2] ;
			n += 3 ;
		}
	}
	//printf("n %d\n",n) ;
	n = 0 ;
	for(j = 0;j < ssy;j++){
		for(i = 0;i < ssx;i++){
			actmap.at<cv::Vec3b>(j, i)[2] = 255 ;
			actmap.at<cv::Vec3b>(j, i)[1] = 255 ;
			actmap.at<cv::Vec3b>(j, i)[0] = 255 ;
		}
	}

	std::chrono::system_clock::time_point start, end ;
	start = std::chrono::system_clock::now() ;

	for( j = sch_y1;j < sch_y2;j += skips_n){
		for(i = sch_x1;i < sch_x2;i += skips_n){
			tgt_x1 = (i < 0)?0:i ;
			tgt_y1 = (j < 0)?0:j ;
			tgt_x2 = (i + isx > ssx)?ssx:i + isx ;
			tgt_y2 = (j + isy > ssy)?ssy:j + isy ;
			tgt_s = (tgt_x2 - tgt_x1) * (tgt_y2 - tgt_y1) ;

			if(count++){

				com_x1 = fmax(ref_x1, tgt_x1) ;
				com_y1 = fmax(ref_y1, tgt_y1) ;
				com_x2 = fmin(ref_x2, tgt_x2) ;
				com_y2 = fmin(ref_y2, tgt_y2) ;
				
				if(com_x2 - com_x1 > 0 && com_y2 - com_y1 > 0){
					com_s = (com_x2 - com_x1) * (com_y2 - com_y1) ;
				}else com_s = 0 ;

				// Non-overlappint parts
				dif_s = tgt_s - com_s ;

				// Upper value
				uplimit = (upl_mag * fmin(com_s, ref_r * ref_s) + dif_s) / (double)tgt_s ;

				if(uplimit < ref_r){

					if(ip_safe(i, j, ssx, ssy)){
						actmap.at<cv::Vec3b>(j, i)[2] = 255 ;
						actmap.at<cv::Vec3b>(j, i)[1] = 255 ;
						actmap.at<cv::Vec3b>(j, i)[0] = 255 ;
					}
				}else{

					// Similarity
					nn++ ;
					rate = calc_match_hist(&SS_color, &II_color, j, i, skips_n, isrc, itmp, t_leng, s_leng) ;
					if(ip_safe(i, j, ssx, ssy)){
						actmap.at<cv::Vec3b>(j, i)[2] = 255 ;
						actmap.at<cv::Vec3b>(j, i)[1] = 0 ;
						actmap.at<cv::Vec3b>(j, i)[0] = 0 ;
					}
					//actmap.at<cv::Vec3b>(tgt_y1, tgt_x1)[2] = 255 ;
					//actmap.at<cv::Vec3b>(tgt_y1, tgt_x1)[1] = 0 ;
					//actmap.at<cv::Vec3b>(tgt_y1, tgt_x1)[0] = 0 ;

					ref_x1 = tgt_x1 ;
					ref_x2 = tgt_x2 ;
					ref_y1 = tgt_y1 ;
					ref_y2 = tgt_y2 ;
					ref_s = tgt_s ;
					ref_r = rate ;

					// Registration of template position
					if(rate > 0.6){
						if(s_rate < ref_r){
							//si = rint(i / image_mag) ;
							//sj = rint(j / image_mag) ;
							si = i ;
							sj = j ;
							s_rate = ref_r ;
							found++ ;
							//printf("s_rate : %lf \n",s_rate) ;
						}
					}
					if(found > fmaxx) break ; // 1024
				}
			}else{
				ref_x1 = tgt_x1 ;
				ref_x2 = tgt_x2 ;
				ref_y1 = tgt_y1 ;
				ref_y2 = tgt_y2 ;
				ref_s = tgt_s ;
				ref_r = calc_match_hist(&SS_color, &II_color, j, i, skips_n, isrc, itmp, t_leng, s_leng) ;
			}
		}
	}

	I_color = cv::imread(ifname, cv::IMREAD_GRAYSCALE);
	S_color = cv::imread(iifname, cv::IMREAD_GRAYSCALE);
	NcorrBy1D( &I_color, &S_color, si, sj ) ;

	end = std::chrono::system_clock::now() ;
	const double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() ;

	for( j=0 ; j<ssy ; j++ ) {
		for( i=0 ; i<ssx ; i++ ) {
			RESR.at<unsigned char>(j,i) = SS_color.at<cv::Vec3b>(j,i)[2] ;
			RESG.at<unsigned char>(j,i) = SS_color.at<cv::Vec3b>(j,i)[1] ;
			RESB.at<unsigned char>(j,i) = SS_color.at<cv::Vec3b>(j,i)[0] ;
		}
	}

	ei = si + isx ;
	ej = sj + isy ;

	cv::rectangle(RESR,cv::Point(si-1,sj-1),cv::Point(ei+1,ej+1),cv::Scalar(255),1);
	cv::rectangle(RESG,cv::Point(si-1,sj-1),cv::Point(ei+1,ej+1),cv::Scalar(0),1);
	cv::rectangle(RESB,cv::Point(si-1,sj-1),cv::Point(ei+1,ej+1),cv::Scalar(0),1);

	cv::rectangle(RESR,cv::Point(si,sj),cv::Point(ei,ej),cv::Scalar(200),1);
	cv::rectangle(RESG,cv::Point(si,sj),cv::Point(ei,ej),cv::Scalar(0),1);
	cv::rectangle(RESB,cv::Point(si,sj),cv::Point(ei,ej),cv::Scalar(255),1);

	cv::rectangle(RESR,cv::Point(si+1,sj+1),cv::Point(ei-1,ej-1),cv::Scalar(255),1);
	cv::rectangle(RESG,cv::Point(si+1,sj+1),cv::Point(ei-1,ej-1),cv::Scalar(0),1);
	cv::rectangle(RESB,cv::Point(si+1,sj+1),cv::Point(ei-1,ej-1),cv::Scalar(0),1);

	//cv::putText(RESR,std::to_string(n),cv::Point(si,sj-4),cv::FONT_HERSHEY_SIMPLEX,0.4,cv::Scalar(255));
	//cv::putText(RESG,std::to_string(n),cv::Point(si,sj-4),cv::FONT_HERSHEY_SIMPLEX,0.4,cv::Scalar(255));
	//cv::putText(RESB,std::to_string(n),cv::Point(si,sj-4),cv::FONT_HERSHEY_SIMPLEX,0.4,cv::Scalar(0));

	std::vector<cv::Mat> color_img;
	cv::Mat color_img_dst;
	color_img.push_back(RESB);
	color_img.push_back(RESG);
	color_img.push_back(RESR);
	cv::merge(color_img,color_img_dst);
	cv::imwrite(resfname,color_img_dst);

	fp = fopen( ssfname, "w") ;
	fprintf(fp, "time %lf[ms]\n", time) ;
	fprintf(fp, "%d\n", 1) ;
	fprintf(fp, "score %d\n", 0) ;
	fprintf(fp, "i %3d, j %3d\n", si, sj ) ;
	fclose( fp ) ;
	cv::imwrite(mapfname, actmap) ;
	
	II_color.release();
	SS_color.release();
	RESR.release();
	RESG.release();
	RESB.release();
	//actmap.release() ;
	free(isrc) ;
	free(itmp) ;

	//printf("found number : %d\n",found) ;
	//printf("match number : %d\n",nn) ;

	return 0 ;

}

double calc_match_rate(cv::Mat *SS_color, cv::Mat *II_color, int j, int i, int skip_n, unsigned char *isrc, unsigned char *itmp, int t_leng, int s_leng){

	double error, re, ge, be, th ;
	int isx, isy ,ssx, ssy ;
	int sch_x1, sch_y1, sch_x2, sch_y2 ;
	int match, pixel, jj, ii ; 

	isx = II_color->cols ;
	isy = II_color->rows ;
	ssx = SS_color->cols ;
	ssy = SS_color->rows ;
	sch_x1 = (i < 0)?abs(i):0 ;
	sch_y1 = (j < 0)?abs(j):0 ;
	sch_x2 = (i + isx > ssx)?(ssx - i):isx ;
	sch_y2 = (j + isy > ssy)?(ssy - j):isy ;

	match = 0 ;
	pixel = 0 ;
	th = 14.0 ;

	for(jj = sch_y1;jj < sch_y2;jj += skip_n){
		for(ii = sch_x1;ii < sch_x2;ii += skip_n){
			if(ip_safe(i + ii, j + jj, ssx, ssy)){
				pixel++ ;
				re = abs(isrc[XY1D((i + ii) * 3, j + jj, s_leng) + 2] - itmp[XY1D(ii * 3, jj, t_leng) + 2]) ;
				ge = abs(isrc[XY1D((i + ii) * 3, j + jj, s_leng) + 1] - itmp[XY1D(ii * 3, jj, t_leng) + 1]) ;
				be = abs(isrc[XY1D((i + ii) * 3, j + jj, s_leng) + 0] - itmp[XY1D(ii * 3, jj, t_leng) + 0]) ;
				error = sqrt(re * re + ge * ge + be * be) ;

				if(error < th) match++ ;
			}
		}
	}

	return (double)match / (double)pixel ;

}

double calc_match_hist(cv::Mat *SS_color, cv::Mat *II_color, int j, int i, int skip_n, unsigned char *isrc, unsigned char *itmp, int t_leng, int s_leng){

	double rate, re, ge, be ;
	int isx, isy ,ssx, ssy, th ;
	int sch_x1, sch_y1, sch_x2, sch_y2 ;
	int match, pixel, jj, ii, n ;
	int max_col, rbit ;
	double srhist[512], sghist[512], sbhist[512] ;
	double trhist[512], tghist[512], tbhist[512] ;

	isx = II_color->cols ;
	isy = II_color->rows ;
	ssx = SS_color->cols ;
	ssy = SS_color->rows ;
	sch_x1 = (i < 0)?abs(i):0 ;
	sch_y1 = (j < 0)?abs(j):0 ;
	sch_x2 = (i + isx > ssx)?(ssx - i):isx ;
	sch_y2 = (j + isy > ssy)?(ssy - j):isy ;

	match = 0 ;
	pixel = 0 ;
	th = 1 ; // Allowed search range for color or gray
	rbit = 0 ; // max 8
	max_col	= pow(2, 8 - rbit) ;

	for(n = 0;n < 512;n++){
		srhist[n] = sghist[n] = sbhist[n] = 0.0 ;
		trhist[n] = tghist[n] = tbhist[n] = 0.0 ;
	}

	// Quantize 256 steps by number of bits (0-256)
	// Make histograms R,G,B
	for(jj = sch_y1;jj < sch_y2;jj += skip_n){
		for(ii = sch_x1;ii < sch_x2;ii += skip_n){
			if(ip_safe(i + ii, j + jj, ssx, ssy)){

				pixel++ ;
				srhist[isrc[XY1D((i+ii)*3,j+jj,s_leng)+2]>>rbit]+=1;
				sghist[isrc[XY1D((i+ii)*3,j+jj,s_leng)+1]>>rbit]+=1;
				sbhist[isrc[XY1D((i+ii)*3,j+jj,s_leng)+0]>>rbit]+=1;
				trhist[itmp[XY1D(ii*3,jj,t_leng)+2]>>rbit]+=1;
				tghist[itmp[XY1D(ii*3,jj,t_leng)+1]>>rbit]+=1;
				tbhist[itmp[XY1D(ii*3,jj,t_leng)+0]>>rbit]+=1;

				//printf("R %d \n",isrc[XY1D((i+ii)*3,j+jj,s_leng)+2]>>rbit);
				//printf("G %d \n",isrc[XY1D((i+ii)*3,j+jj,s_leng)+1]);
				//printf("B %d \n",isrc[XY1D((i+ii)*3,j+jj,s_leng)+0]);
			}
		}
	}

	// Normalization
	rate = 0.0 ;
	for(n = 0;n < max_col;n++){
		//printf("%lf %lf %lf %d \n",srhist[n],sghist[n],sbhist[n],n) ;
		//printf("%lf %lf %lf %d \n",trhist[n],tghist[n],tbhist[n],n) ;
		srhist[n] = srhist[n] / (double)pixel ;
		sghist[n] = sghist[n] / (double)pixel ;
		sbhist[n] = sbhist[n] / (double)pixel ;
		trhist[n] = trhist[n] / (double)pixel ;
		tghist[n] = tghist[n] / (double)pixel ;
		tbhist[n] = tbhist[n] / (double)pixel ;
	}

	// Histogram intersection
	for(n = 0;n < max_col;n++){
		re = fmin(srhist[n], trhist[n]) ;
		ge = fmin(sghist[n], tghist[n]) ;
		be = fmin(sbhist[n], tbhist[n]) ;
		rate  += (re + ge + be) / 3.0 ;
		//printf("rate %lf n %d pixel %d \n",rate, n, pixel) ;
		//printf("re %lf ge %lf be %lf\n",re,ge,be) ;
	}

	return rate ;

}

void NcorrBy1D( cv::Mat *I, cv::Mat *S, int si, int sj ){

	register int	i, j, ii, jj ;
	int				isx, isy, ssx, ssy, ssi, ssj ;
	int				Ndata, sfg, sff, sf, sgg, sg, area ;
	int				f, g ;
	double			enume, darea, corr ;
	double			divF, divT, n ;


	isx = I->cols;
	isy = I->rows;
	ssx = S->cols ;
	ssy = S->rows ;
	n = -1 ;

	printf("si : %d, sj : %d\n",si, sj) ;


	// Normalizing template
	area = isx * isy  ;
	darea = (double)area ;
	sgg = sg = 0 ;
	for( jj=0 ; jj<isy ; jj++ ) {
		for( ii=0 ; ii<isx ; ii++ ) {
			g    = I->at<unsigned char>(jj,ii) ;
			sgg += g * g ;
			sg	+= g ;
		}
	}
	divT = darea * (double)sgg - (double)sg * (double)sg ;

	// Generate score map
	sf = sff = sfg = 0 ;
	for( jj=sj - 40 ; jj<sj + 50 ; jj++ ) {
		for( ii=si - 40 ; ii<si + 50 ; ii++ ) {
			sfg = sff = sf = 0 ;
			for( j=0 ; j<isy ; j++ ) {
				for( i=0 ; i<isx ; i++ ) {
					f = (int)( S->at<unsigned char>(jj+j,ii+i) ) ;
					sf	+= f ;
					sff += f * f ;
					sfg += f * (int)(I->at<unsigned char>(j,i)) ;
				}
			}
			divF  = darea * (double)sff - (double)sf * (double)sf ;
			enume = darea * (double)sfg - (double)sf * (double)sg ;
			corr  = enume / sqrt( divF * divT ) ;
			if( divF == 0.0 || divT == 0.0) corr = 0.0;

			if(n < corr){
				ssi = ii ;
				ssj = jj ;
				n = corr ;
				//printf("corr : %lf\n",n) ;
			}
		}
	}

	si = ssi ;
	sj = ssj ;
	printf("ssi : %d, ssj : %d\n",si, sj) ;

	return ;

}