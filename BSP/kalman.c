#include "kalman.h"
#include "mpu6050.h"
#include "inv_mpu.h"
#include "inv_mpu_dmp_motion_driver.h"
//���������㷨��

short aacx,aacy,aacz;		//���ٶȴ�����ԭʼ���� 
short gyrox,gyroy,gyroz;	//������ԭʼ���� 
short temperature;			//�������¶�����
float accel_x;	     		//X����ٶ�ֵ
float accel_y;	    		//Y����ٶ�ֵ
float accel_z;	     		//Z����ٶ�ֵ
float gyro_x;						//X�������ǽ��ٶ�����
float gyro_y;        		//Y�������ǽ��ٶ�����
float gyro_z;		 				//Z�������ǽ��ٶ�����
float yaw_raw;  		 		//���߽�yawԭʼ����
float yaw_kalman; 	 		//yaw�˲�������
float pitch_raw;  		 	//������pitchԭʼ����
float pitch_kalman; 	 	//pitch�˲�������
float roll_raw;  		   	//�����rollԭʼ����
float roll_kalman; 		 	//roll�˲�������					

int16_t temp;        //�¶�
int16_t gx,gy,gz;    //������ٶ�
int16_t ax,ay,az;    //������ٶ�
//�ж����ݽ����Ƿ�����
int8_t temp_ret;
int8_t a_ret;
int8_t g_ret;

//���������ݼ���
void Angle_Cal(void)
{
	//1. ԭʼ���ݻ�ȡ
//	float accx,accy,accz,acctemp;//������Ǽ��ٶ�ֵ
	//��ȡ�¶�
	temp_ret=MPU_Get_Temperature();						//�õ��¶���Ϣ
	//��ȡ���ٶȴ���������			
  a_ret =  MPU_Get_Accelerometer(&ax,&ay, &az);	;
  //������������� 
  g_ret=MPU_Get_Gyroscope(&gx, &gy, &gz);		//�õ�����������
	accel_x = ax;//x����ٶ�ֵ�ݴ�
	accel_y = ay;//y����ٶ�ֵ�ݴ�
	accel_z = az;//z����ٶ�ֵ�ݴ�
	gyro_x  = gx;//x��������ֵ�ݴ�
	gyro_y  = gy;//y��������ֵ�ݴ�
	gyro_z  = gz;//z��������ֵ�ݴ�
	
//	//2.�Ǽ��ٶ�ԭʼֵ�������	
//	//���ٶȴ��������üĴ���0X1C��д��0x01,���÷�ΧΪ��2g�������ϵ��2^16/4 = 16384LSB/g
//	if(accel_x<32764) accx=accel_x/16384.0;//����x����ٶ�
//	else              accx=1-(accel_x-49152)/16384.0;
//	if(accel_y<32764) accy=accel_y/16384.0;//����y����ٶ�
//	else              accy=1-(accel_y-49152)/16384.0;
//	if(accel_z<32764) accz=accel_z/16384.0;//����z����ٶ�
//	else              accz=(accel_z-49152)/16384.0;
//	//���ٶȷ����й�ʽ�����������ˮƽ������ϵ֮��ļн�
//	acctemp=sqrt(accx*accx+accy*accy);
//	pitch_raw=(atan(accy/accz))*180/3.14;
//	roll_raw=(atan(accx/accz))*180/3.14;
//  yaw_raw=(atan(acctemp/accz))*180/3.14;
//	//�жϼ����Ƕȵ�������											
//	if(accel_x<32764) roll_raw = +roll_raw;
//	if(accel_x>32764) roll_raw = -roll_raw;
//	if(accel_y<32764) pitch_raw = +pitch_raw;
//	if(accel_y>32764) pitch_raw = -pitch_raw;
//	if(accel_z<32764) yaw_raw = +yaw_raw;
//	if(accel_z>32764) yaw_raw = -yaw_raw;
	while(mpu_dmp_get_data(&pitch_raw, &roll_raw, &yaw_raw));	
	//3.���ٶ�ԭʼֵ�������
	//���������üĴ���0X1B��д��0x18�����÷�ΧΪ��2000deg/s�������ϵ��2^16/4000=16.4LSB/(��/S)
	////������ٶ�
	if(gyro_x<32768) gyro_x=-(gyro_x/16.4);
	if(gyro_x>32768) gyro_x=+(65535-gyro_x)/16.4;
	if(gyro_y<32768) gyro_y=-(gyro_y/16.4);
	if(gyro_y>32768) gyro_y=+(65535-gyro_y)/16.4;
	if(gyro_z<32768) gyro_z=-(gyro_z/16.4);
	if(gyro_z>32768) gyro_z=+(65535-gyro_z)/16.4;

	//4.���ÿ���������
	Kalman_Cal_Pitch(pitch_raw,gyro_x);  //�������˲�����pitch
	Kalman_Cal_Roll(roll_raw,gyro_y);  //�������˲�����roll
	Kalman_Cal_Yaw(yaw_raw,gyro_z); //�������˲�����yaw


}
 
//����������		
static float Q_angle = 0.001;		//�Ƕ��������Ŷȣ��Ƕ�������Э����
static float Q_gyro  = 0.003;		//���ٶ��������Ŷȣ����ٶ�������Э����  
static float R_angle = 0.5;		//���ٶȼƲ���������Э����
static float dt      = 0.01;		//�������ڼ�������������10ms

void Kalman_Cal_Pitch(float acc,float gyro) //�������˲�pitch�����		
{
	static float Q_bias;	//Q_bias:�����ǵ�ƫ��
	static float K_0, K_1;	//����������  K_0:���ڼ������Ź���ֵ  K_1:���ڼ������Ź���ֵ��ƫ�� t_0/1:�м����
	static float PP[2][2] = { { 1, 0 },{ 0, 1 } };//����Э�������P����ʼֵΪ��λ��

	/*
		�������˲���ʹ�ò���
		(1) ѡ��״̬�����۲���
		(2) ��������
		(3) ��ʼ������
		(4) ���빫ʽ����
		(5) ���ڳ�����P��Q
	*/
	/*
	X(k)��kʱ��ϵͳ״̬				Z(k)��kʱ�̲���ֵ
	U(k)��kʱ�̶�ϵͳ������		H������ϵͳ����
																					 ����
	A/F��״̬ת�ƾ���					W(k)���������� ----> Q
																					 ����	
	B�����ƾ���								V(k)���������� ----> R
	
										��ɢ����ϵͳ
	ϵͳ������X(k|k-1) = AX(k-1|k-1) + BU(k) + (W(k))
	����ֵ��Z(k) = HX(k) + V(k)
	*/
	/*
	1. �������
* * *��ʽ1��X(k|k-1) = AX(k-1|k-1) + BU(k) + (W(k))

		X = (Angle,Q_bias)
		A(1,1) = 1,A(1,2) = -dt
		A(2,1) = 0,A(2,2) = 1
		ע����������[���������
		Ԥ�⵱ǰ�Ƕ�ֵ��
		[ angle ] 	[1 -dt][ angle ]	 [dt]
		[ Q_bias] = [0  1 ][ Q_bias] + [ 0] * newGyro(���ٶȼƲ���ֵ)
		��
		angle = angle - Q_bias*dt + newGyro * dt
		Q_bias = Q_bias
	*/
	pitch_kalman += (gyro - Q_bias) * dt; //״̬����,�Ƕ�ֵ�����ϴ����ŽǶȼӽ��ٶȼ���Ư�����
	
	/*
	2. Ԥ��Э�������
* * *��ʽ2��P(k|k-1)=AP(k-1|k-1)A^T + Q 

		�����������ϵͳ����
				[1 -dt]
		A = [0  1 ]
		
		ϵͳ����Э��������Q���壺
		| cov(angle,angle)  cov(Q_bias,angle) |
		|	cov(angle,Q_bias) cov(Q_bias,Q_bias)|
			 �Ƕ������ͽ��ٶ�Ư�������໥����
		| D( angle )  			0 	 |
	= |			0 			D( Q_bias )|
		��Q_angle��Q_bias�ķ���Ϊ������
		���ɾ�������ó�
		
		��D( angle )  = Q_angle 
			D( Q_bias ) = Q_gyro 
			
		����һ��Ԥ��Э�������ΪP(k-1)
											|a(k-1)  b(k-1)|
											|c(k-1)  d(k-1)|
		����Ԥ��Э�������P(k)
											|a(k)  b(k)|
		                  |c(k)  d(k)|
		�ɹ�ʽ2��P(k|k-1)=AP(k-1|k-1)A^T + Q ��
		|a(k)  b(k)|		|1 -dt|	|a(k-1) b(k-1)| |1 	 0|		| D( angle )  			0 	 |
		|c(k)  d(k)| =  |0  1 | |c(k-1) d(k-1)| |-dt 1| + |			0 			D( Q_bias )|
	
		��һ����
		|a(k)  b(k)|		|a(k-1) - [c(k-1) + b(k-1)]*dt + d(dt)^2		b(k-1) - d(k-1)*dt|		| D( angle )  			0 	 |
	  |c(k)  d(k)| =  |				c(k-1) - d(k-1)*dt													d(k-1)		| + |			0 			D( Q_bias )|
		
		����dt^2̫С����dt^2ʡ��
	*/
	
	PP[0][0] = PP[0][0] + Q_angle - (PP[0][1] + PP[1][0])*dt;
	PP[0][1] = PP[0][1] - PP[1][1]*dt;
	PP[1][0] = PP[1][0] - PP[1][1]*dt;
	PP[1][1] = PP[1][1] + Q_gyro;
	
	/*
		3. ������������
			ϵͳ�������� Z(k) = HX(k) + V(k)
			ϵͳ����ϵ�� H = [1, 0]
			��Ϊ����������Դ�����
			����
			measure = newAngle
	*/
	
	/*
		4. ���㿨��������
* * *��ʽ3��Kg(k)= P(k|k-1)H^T/(HP(k|k-1)H^T+R)
				Kg = (K_0,K_1) ��Ӧangle,Q_bias����
				H = (1,0)
	*/
	K_0 = PP[0][0] / (PP[0][0] + R_angle);
	K_1 = PP[1][0] / (PP[0][0] + R_angle);

	/*
		5. ���㵱ǰ���Ż�����ֵ
* * *��ʽ4��X(k|k) = X(k|k-1) + kg(k)[z(k) - HX(k|k-1)]
		angle = angle + K_0*(newAngle - angle)
		Q_bias = Q_bias + K_1*(newAngle - angle)
	*/
		
	pitch_kalman = pitch_kalman + K_0 * (acc - pitch_kalman);
	Q_bias = Q_bias + K_1 * (acc - pitch_kalman);

	/*
		6. ����Э�������
* * *��ʽ5��P(k|k)=[I-Kg(k)H]P(k|k-1)
	*/
	PP[0][0] = PP[0][0] - K_0 * PP[0][0];
	PP[0][1] = PP[0][1] - K_0 * PP[0][1];
	PP[1][0] = PP[1][0] - K_1 * PP[0][0];
	PP[1][1] = PP[1][1] - K_1 * PP[0][1];

}

void Kalman_Cal_Roll(float acc,float gyro) //�������˲�roll�����				
{
	static float Q_bias;	//Q_bias:�����ǵ�ƫ��  Angle_err:�Ƕ�ƫ�� 
	static float K_0, K_1;	//����������  K_0:���ڼ������Ź���ֵ  K_1:���ڼ������Ź���ֵ��ƫ�� t_0/1:�м����
	static float PP[2][2] = { { 1, 0 },{ 0, 1 } };//����Э�������P����ʼֵΪ��λ��
	roll_kalman += (gyro - Q_bias) * dt; //״̬����,�Ƕ�ֵ�����ϴ����ŽǶȼӽ��ٶȼ���Ư�����
	PP[0][0] = PP[0][0] + Q_angle - (PP[0][1] + PP[1][0])*dt;
	PP[0][1] = PP[0][1] - PP[1][1]*dt;
	PP[1][0] = PP[1][0] - PP[1][1]*dt;
	PP[1][1] = PP[1][1] + Q_gyro;
	K_0 = PP[0][0] / (PP[0][0] + R_angle);
	K_1 = PP[1][0] / (PP[0][0] + R_angle);
	roll_kalman = roll_kalman + K_0 * (acc - roll_kalman);
	Q_bias = Q_bias + K_1 * (acc - roll_kalman);
	PP[0][0] = PP[0][0] - K_0 * PP[0][0];
	PP[0][1] = PP[0][1] - K_0 * PP[0][1];
	PP[1][0] = PP[1][0] - K_1 * PP[0][0];
	PP[1][1] = PP[1][1] - K_1 * PP[0][1];
}

void Kalman_Cal_Yaw(float acc,float gyro) //�������˲�yaw�����				
{
	static float Q_bias;	//Q_bias:�����ǵ�ƫ��  Angle_err:�Ƕ�ƫ�� 
	static float K_0, K_1;	//����������  K_0:���ڼ������Ź���ֵ  K_1:���ڼ������Ź���ֵ��ƫ�� t_0/1:�м����
	static float PP[2][2] = { { 1, 0 },{ 0, 1 } };//����Э�������P����ʼֵΪ��λ��
	roll_kalman += (gyro - Q_bias) * dt; //״̬����,�Ƕ�ֵ�����ϴ����ŽǶȼӽ��ٶȼ���Ư�����
	PP[0][0] = PP[0][0] + Q_angle - (PP[0][1] + PP[1][0])*dt;
	PP[0][1] = PP[0][1] - PP[1][1]*dt;
	PP[1][0] = PP[1][0] - PP[1][1]*dt;
	PP[1][1] = PP[1][1] + Q_gyro;
	K_0 = PP[0][0] / (PP[0][0] + R_angle);
	K_1 = PP[1][0] / (PP[0][0] + R_angle);
	yaw_kalman = yaw_kalman + K_0 * (acc - yaw_kalman);
	Q_bias = Q_bias + K_1 * (acc - yaw_kalman);
	PP[0][0] = PP[0][0] - K_0 * PP[0][0];
	PP[0][1] = PP[0][1] - K_0 * PP[0][1];
	PP[1][0] = PP[1][0] - K_1 * PP[0][0];
	PP[1][1] = PP[1][1] - K_1 * PP[0][1];
}
