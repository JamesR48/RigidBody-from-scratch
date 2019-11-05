#ifdef _WIN32
#include "glut.h"
#elif defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
#include <Eigen/Dense>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <vector>
using Eigen::Vector4f;
using Eigen::Matrix4f;
using Eigen::MatrixXd;

using namespace std;
using namespace Eigen;
#define PI 3.14159265

bool M = GL_TRUE; //Malla (Tecla Z)
bool E = GL_TRUE; //Eje (Tecla X)
bool Q = GL_TRUE; //Figura (Tecla Q)
bool N = GL_TRUE; //Normales (Tecla N)
bool V = GL_TRUE; //Vertices (Tecla V)
bool C = GL_TRUE; //Caras (Tecla C)
bool B = GL_TRUE; //Aristas (Tecla B)

int cuentaAuilixar = 0;

GLfloat cop_x;
GLfloat cop_y;
GLfloat cop_z;
GLfloat pref_x;
GLfloat pref_y;
GLfloat pref_z;
GLfloat znear = 1.0;
GLfloat zfar = 300;
GLfloat xWmax = 500;
GLfloat yWmax = 500;

Vector3f pos;                     //Posicion cm
Vector3f vel;                     //Velocidad
Vector3f F;                       //Fuerza
Vector3f P´[12];                  //Posiciones particulas
Vector3f Particula[12];           //Posiciones globales particulas
Vector4f Pref, Punto;
Matrix3f I_0;                     //Tensor de inercia
Matrix3f I_inv;                   //Tensor de inercia inversa
Vector3f CM, CM´;                 //Centro de masas (coord. locales , coord. globales)
Vector3f T;                       //Torque
Vector3f L;                       //Momento lineal
Vector3f W;                       //Velocidad angular
Matrix3f W´;                      //Aceleracion angular
Matrix3f Rot, R´;                 //Matriz de rotacion & Derivada de la rotacion

Matrix4f Rx, Ry, Rz, Traslacion;

GLint anchoG = 100, pasoG = 2;
GLint nvertices = 6;
bool CO = false;

float m = 1;      //Masa
float h = 0.0066;   //Salto

#define PI 3.14159265

float alfa = -45;
float beta = 20;
float gama = 30;
float x;               //Coordenadas de posicion centro de masa
float y;
float z;
float vx = 0;          //Coordenadas de la velocidad centro de masa 
float vy = -3.7;
float vz = 0;
float X1 = 4;
float Y1 = 5;
float Z1 = 0;

class Listavertices {
private:
	Vector3f vector;
	float verticex, verticey, verticez;
public:
	double x();
	double y();
	double z();
	void setx(float x_);
	void sety(float);
	void setz(float);
	Vector3f retornarvector();
	Listavertices();//CONSTRUCTOR
	Listavertices(float x2, float y2, float z2);
};

Listavertices::Listavertices() {

}

Listavertices::Listavertices(float x2, float y2, float z2) {
	verticex = x2;
	verticey = y2;
	verticez = z2;
}

double Listavertices::x() {
	return verticex;
}

double Listavertices::y() {
	return verticey;
}

double Listavertices::z() {
	return verticez;
}

void Listavertices::setx(float x_)
{
	verticex = x_;
}

void Listavertices::sety(float y_)
{
	verticey = y_;
}

void Listavertices::setz(float z_)
{
	verticez = z_;
}

Vector3f Listavertices::retornarvector()
{
	vector << verticex, verticey, verticez;
	return vector;
}


Listavertices listaP[12];


void InicializarPoligono()
{
	listaP[0] = Listavertices(2, 0, 0);
	listaP[1] = Listavertices(1, 0, -1.73205);
	listaP[2] = Listavertices(-1, 0, -1.73205);
	listaP[3] = Listavertices(-2, 0, 0);
	listaP[4] = Listavertices(-1, 0, 1.73205);
	listaP[5] = Listavertices(1, 0, 1.73205);
	listaP[6] = Listavertices(2, 4, 0);
	listaP[7] = Listavertices(1, 4, -1.73205);
	listaP[8] = Listavertices(-1, 4, -1.73205);
	listaP[9] = Listavertices(-2, 4, 0);
	listaP[10] = Listavertices(-1, 4, 1.73205);
	listaP[11] = Listavertices(1, 4, 1.73205);
}


void CentrodeMasa()
{
	InicializarPoligono();
	CM = (listaP[0].retornarvector() + listaP[1].retornarvector() + listaP[2].retornarvector() + listaP[3].retornarvector() +
		listaP[4].retornarvector() + listaP[5].retornarvector() + listaP[6].retornarvector() + listaP[7].retornarvector() +
		listaP[8].retornarvector() + listaP[9].retornarvector() + listaP[10].retornarvector() + listaP[11].retornarvector()) / 12;
	P´[0] = listaP[0].retornarvector() - CM;
	P´[1] = listaP[1].retornarvector() - CM;
	P´[2] = listaP[2].retornarvector() - CM;
	P´[3] = listaP[3].retornarvector() - CM;
	P´[4] = listaP[4].retornarvector() - CM;
	P´[5] = listaP[5].retornarvector() - CM;
	P´[6] = listaP[6].retornarvector() - CM;
	P´[7] = listaP[7].retornarvector() - CM;
	P´[8] = listaP[8].retornarvector() - CM;
	P´[9] = listaP[9].retornarvector() - CM;
	P´[10] = listaP[10].retornarvector() - CM;
	P´[11] = listaP[11].retornarvector() - CM;
	CM´ = (P´[0] + P´[1] + P´[2] + P´[3] + P´[4] + P´[5] + P´[6] + P´[7] + P´[8] + P´[9] + P´[10] + P´[11]) / 12;
	x = CM´.x();
	y = CM´.y();
	z = CM´.z();
	pref_x = pos.x();
	pref_y = pos.y();
	pref_z = pos.z();
	cop_x = pos.x();
	cop_y = pos.y() + 5;
	cop_z = pos.z() + 50;
	//cout << "Centro de masa local:" << endl << CM << endl;
	//cout << "Centro de masa global:" << endl << CM´ << endl;
}

void Tensordeinercia() {
	Matrix3f I[12];
	for (int i = 0; i < 12; i++)
	{

		I[i] << m * ((pow(P´[i].y(), 2)) + (pow(P´[i].z(), 2))), -m * P´[i].x()*P´[i].y(), -m * P´[i].x()*P´[i].z(),

			-m * P´[i].y()*P´[i].x(), m*((pow(P´[i].x(), 2)) + (pow(P´[i].z(), 2))), -m * P´[i].y()*P´[i].z(),

			-m * P´[i].z()*P´[i].x(), -m * P´[i].z()*P´[i].y(), m*((pow(P´[i].x(), 2)) + (pow(P´[i].y(), 2)));
	}

	for (int j = 0; j < 12; j++)
	{
		I_0 += I[j];
	}
	//cout << "Momento de inercia; " << endl << I_0 << endl;
	I_inv = I_0.inverse();

}

void DibujarObjeto()
{
	InicializarPoligono();

	/*for (int i = 0; i < 12; i++){

		Particula[i] = (Rot*listaP[i].retornarvector() + CM);

	}


	for (int k = 0; k <= nvertices; k++)

	{

		for (int j = 0; j < nvertices - 1; j++){


			glBegin(GL_LINE_LOOP);

			glColor3f(1, 1, 1);

			glVertex3f(Particula[j].x(), Particula[j].y(), Particula[j].z());

			glVertex3f(Particula[j + 1].x(), Particula[j + 1].y(), Particula[j + 1].z());

			glVertex3f(Particula[j + 7].x(), Particula[j + 7].y(), Particula[j + 7].z());

			glVertex3f(Particula[j + 6].x(), Particula[j + 6].y(), Particula[j + 6].z());

			glEnd();

		}


		glBegin(GL_LINE_LOOP);

		glColor3f(1, 1, 1);

		glVertex3f(Particula[5].x(), Particula[5].y(), Particula[5].z());

		glVertex3f(Particula[0].x(), Particula[0].y(), Particula[0].z());

		glVertex3f(Particula[6].x(), Particula[6].y(), Particula[6].z());

		glVertex3f(Particula[11].x(), Particula[11].y(), Particula[11].z());

		glEnd();

	}


	for (int i = 0; i < 12; i++){

		Particula[i] = (Rot*listaP[i].retornarvector() + CM);

	}*/

	for (int k = 0; k <= nvertices; k++)
	{
		for (int j = 0; j < nvertices - 1; j++) {

			glBegin(GL_LINE_LOOP);

			glColor3f(1, 1, 1);

			glVertex3f(P´[j].x(), P´[j].y(), P´[j].z());

			glVertex3f(P´[j + 1].x(), P´[j + 1].y(), P´[j + 1].z());

			glVertex3f(P´[j + 7].x(), P´[j + 7].y(), P´[j + 7].z());

			glVertex3f(P´[j + 6].x(), P´[j + 6].y(), P´[j + 6].z());

			glEnd();

		}


		glBegin(GL_LINE_LOOP);

		glColor3f(1, 1, 1);

		glVertex3f(P´[5].x(), P´[5].y(), P´[5].z());

		glVertex3f(P´[0].x(), P´[0].y(), P´[0].z());

		glVertex3f(P´[6].x(), P´[6].y(), P´[6].z());

		glVertex3f(P´[11].x(), P´[11].y(), P´[11].z());

		glEnd();

	}


	glPointSize(5);

	glBegin(GL_POINTS);

	glColor3f(1, 1, 1);

	glVertex3f(CM´.x(), CM´.y(), CM´.z());

	glEnd();


}


void Inicializar(float x0, float y0, float z0, float vx0, float vy0, float vz0)
{
	pos << x0, y0, z0;
	vel << vx0, vy0, vz0;
	L << 0, 0, 0;
	Rot << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;
}

void init(void)
{
	Inicializar(x, y, z, vx, vy, vz);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_FLAT);
	glEnable(GL_DEPTH_TEST);
}

Vector3f GetPos() {
	return pos;
}

Vector3f GetVel() {
	return vel;
}

void limpiar()
{
	F << 0, 0, 0;
}

void AgregarFuerza(Vector3f f0)
{
	F += f0;
}

void Rotacionanalitica()
{
	alfa = alfa * (PI / 180);
	beta = beta * (PI / 180);
	gama = gama * (PI / 180);
	Matrix4f Traslacion, Rc;
	Matrix3f Rx, Ry, Rz, vec;
	Vector3f V, Punto, Punto2, cop, camara;
	Punto << x, y, z;

	Traslacion << 1, 0, 0, -X1,
		0, 1, 0, -Y1,
		0, 0, 1, -Z1,
		0, 0, 0, 1;

	Rx << 1, 0, 0,
		0, cos(alfa), -sin(alfa),
		0, sin(alfa), cos(alfa);

	Ry << cos(beta), 0, sin(beta),
		0, 1, 0,
		-sin(beta), 0, cos(beta);

	Rz << cos(gama), -sin(gama), 0,
		sin(gama), cos(gama), 0,
		0, 0, 1;
	//Punto = Traslacion*Punto;
	//vec = Rz*Ry*Rx;
	Punto = Rz * Ry*Rx*Punto;
	//Punto = Traslacion.inverse()*Punto;
}

void Torque()
{
	Vector3f Ti[12], t;
	Vector3f r[12];
	for (int ro = 0; ro < 12; ro++)
	{
		r[ro] << 0, 0, 2;
	}
	t << 0, 25, 0;
	for (int i = 0; i < 12; i++)
	{
		Ti[i] = (Rot*r[i]).cross(t);
	}
	for (int j = 0; j < 12; j++)
	{
		T += Ti[j];
	}
	cout << "Torque: " << endl << T << endl;
	glPointSize(5);
	glBegin(GL_POINTS);
	glColor3f(1, 1, 1);
	glVertex3f(t.x(), t.y(), t.z());
	glEnd();
}

void Translacion()
{
	Vector3f vnew, posnew;
	Vector3f a = F / m;
	vnew = vel + a * h;
	posnew = pos + vnew * h;
	vel = vnew;
	pos = posnew;
	if (vel[1] < 0 && vel[1]>-2) {
		CO = true;
	}
	//cout << "Posicion: " <<endl << pos << endl;
	//cout << "Velocidad: " <<endl << vel << endl;
}

void Momentoangular()
{
	Vector3f Lnew;
	Lnew = L + h * T;
	L = Lnew;
	//cout << "Momento angular: " <<endl << L << endl;
}

void Velocidadangular()
{
	W = I_inv * L;
	W´ << 0, -W.z(), W.y(),
		W.z(), 0, -W.x(),
		-W.y(), W.x(), 0;
	//cout << "Momento de inercia INV: " <<endl << I_inv << endl;
	//cout << "Velocidad angular: " <<endl << W << endl;
	//cout << "Aceleracion angular: " <<endl << W´ << endl;
}

void Derivadadelarotacion()
{
	Vector3f Rcol1, Rcol2, Rcol3;
	//Rcol1 = W.cross(Rot.col(0));
	//Rcol2 = W.cross(Rot.col(1));
	//Rcol3 = W.cross(Rot.col(2));
	//Rcol1 = W´*Rot.col(0);
	//Rcol2 = W´*Rot.col(1);
	//Rcol3 = W´*Rot.col(2);
	R´ = W´ * Rot;
	/*R´ << Rcol1.x(), Rcol2.x(), Rcol3.x(),
		  Rcol1.y(), Rcol2.y(), Rcol3.y(),
		  Rcol1.z(), Rcol2.z(), Rcol3.z();*/
		  //cout << "Derivada de la rotacion: " <<endl << R´ << endl;
		  //cout << "Rotacion" << endl << Rot << endl;
}

void Integralrotacion()
{
	Matrix3f Rnew;
	Rnew = Rot + h * R´;
	Rot = Rnew;
	//cout << "Rotacion" << Rot << endl;
}

void Repararrotacion()
{
	Vector3f C1, C2, C3;
	C1 = Rot.col(0);
	C2 = Rot.col(1);
	C3 = Rot.col(2);
	C1 = C1 / C1.norm();
	C2 = C2 - (C2.transpose()*C1)*C1;
	C2 = C2 / C2.norm();
	C3 = C1.cross(C2);
	Rot << C1.x(), C2.x(), C3.x(),
		C1.y(), C2.y(), C3.y(),
		C1.z(), C2.z(), C3.z();
	//cout << "Rotacion reparada" <<endl << Rot << endl;
}

void Rotar()
{
	Vector3f Rx = Rot.col(0);
	Vector3f Ry = Rot.col(1);
	Vector3f Rz = Rot.col(2);
	//cout << "Rotacion: " << endl << Rot << endl;
	GLfloat R[4][4];
	GLfloat M[16];
	R[0][0] = Rx.x(); R[0][1] = Ry.x(); R[0][2] = Rz.x(); R[0][3] = 0;
	R[1][0] = Rx.y(); R[1][1] = Ry.y(); R[1][2] = Rz.y(); R[1][3] = 0;
	R[2][0] = Rx.z(); R[2][1] = Ry.z(); R[2][2] = Rz.z(); R[2][3] = 0;
	R[3][0] = 0;      R[3][1] = 0; R[3][2] = 0;      R[3][3] = 1;

	M[0] = R[0][0]; M[4] = R[0][1]; M[8] = R[0][2];  M[12] = R[0][3];
	M[1] = R[1][0]; M[5] = R[1][1]; M[9] = R[1][2];  M[13] = R[1][3];
	M[2] = R[2][0]; M[6] = R[2][1]; M[10] = R[2][2]; M[14] = R[2][3];
	M[3] = R[3][0]; M[7] = R[3][1]; M[11] = R[3][2]; M[15] = R[3][3];
	glMultMatrixf(M);
	DibujarObjeto();
}

void RecalcularI()
{
	I_inv = Rot * I_0.inverse()*Rot.transpose();
}

void Colision(int verX_min, int verX_max, int plano_alt, int verZ_min, int verZ_max)
{
	Vector3f aux2 = GetPos();
	//if (pos[0] > verX_min && pos[0] < verX_max && pos[1] - 30 < plano_alt && pos[2] > verZ_min && pos[2] < verZ_max && CO == true)
	//if (aux2[0] > verX_min && aux2[0] < verX_max && aux2[1]< plano_alt && aux2[2] > verZ_min && aux2[2] < verZ_max && CO == true)
	if (aux2[1] - 2.5 < plano_alt)
	{
		vel[1] += 1.9 * -vel[1];
		F[1] = 0;
		CO = false;
	}
	cout << "vel[1]: " << vel[1] << endl;
	cout << "pos[0]: " << aux2[1] << endl;
	cout << "os[2]: " << aux2[2] << endl;
}

void Rotacionfunciones()
{
	//GLfloat R[4][4];
	//GLfloat Rx[4][4];
	/*Identidad
	R[0][0] = 1; R[0][1] = 0; R[0][2] = 0; R[0][3] = 0;
	R[1][0] = 0; R[1][1] = 1; R[1][2] = 0; R[1][3] = 0;
	R[2][0] = 0; R[2][1] = 0; R[2][2] = 1; R[2][3] = 0;
	R[3][0] = 0; R[3][1] = 0; R[3][2] = 0; R[3][3] = 1;*/
	/*
	//Rotacion en x
	R[0][0] = 1; R[0][1] = 0; R[0][2] = 0; R[0][3] = 0;
	R[1][0] = 0; R[1][1] = cos(alfa); R[1][2] = -sin(alfa); R[1][3] = 0;
	R[2][0] = 0; R[2][1] = sin(alfa); R[2][2] = cos(alfa); R[2][3] = 0;
	R[3][0] = 0; R[3][1] = 0; R[3][2] = 0; R[3][3] = 1;*/
	/*GLfloat M[16];
	M[0] = R[0][0]; M[4] = R[0][1]; M[8] = R[0][2];  M[12] = R[0][3];
	M[1] = R[1][0]; M[5] = R[1][1]; M[9] = R[1][2];  M[13] = R[1][3];
	M[2] = R[2][0]; M[6] = R[2][1]; M[10] = R[2][2]; M[14] = R[2][3];
	M[3] = R[3][0]; M[7] = R[3][1]; M[11] = R[3][2]; M[15] = R[3][3];*/
	//glMultMatrixf(M);
	//glColor3f(1, 0.2, 0.2);
	//DibujarObjeto();
	//glutSolidCube(6);
}

void Malla() {
	glLineWidth(1);
	glBegin(GL_LINES);//grilla
	glColor3f(0, 0.5, 0.5);
	for (int i = -anchoG; i <= anchoG; i += pasoG)
	{
		glVertex3f(-anchoG, 0, i); glVertex3f(anchoG, 0, i);
		glVertex3f(i, 0, -anchoG); glVertex3f(i, 0, anchoG);
	}
	glEnd();
}

void Ejes() {
	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex3f(0, 0, 0); glVertex3f(50, 0, 0);
	glColor3f(0, 1, 0);
	glVertex3f(0, 0, 0); glVertex3f(0, 50, 0);
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, 0); glVertex3f(0, 0, 50);
	glEnd();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glColor3f(1.0, 0.0, 0.0);
	glLoadIdentity();
	//Punto de la camara, centro de vista, vector hacia arriba
	gluLookAt(cop_x, cop_y, cop_z, pref_x, pref_y, pref_z, 0.0, 1.0, 0.0);
	if (M)Malla();
	if (E)Ejes();
	CentrodeMasa();
	Tensordeinercia();
	Torque();
	Momentoangular();
	Velocidadangular();
	Derivadadelarotacion();
	Integralrotacion();
	Repararrotacion();
	Colision(-5, 5, -10, -5, 5);
	glBegin(GL_QUADS);
	glColor3f(1, 0, 0);
	glVertex3f(5, -10, 5);
	glVertex3f(-5, -10, 5);
	glVertex3f(-5, -10, -5);
	glVertex3f(5, -10, -5);
	glEnd();
	Vector3f aux2 = GetPos();
	glTranslatef(aux2[0], aux2[1], aux2[2]);
	Rotar();
	RecalcularI();
	glFlush();
}

void OnTimerGL(int id)
{
	Vector3f g, f1, v, t, t2;
	g << 0, -9.8, 0;
	t << 2, 0, 0;
	float beta = 0.5; //Resistencia del aire
	limpiar();
	AgregarFuerza(g);
	Translacion();
	glutPostRedisplay();
	glutTimerFunc(1, OnTimerGL, 1);
}

void RotacionenY(float theta)
{
	Punto << cop_x, cop_y, cop_z, 1;
	Ry << cos(theta), 0, sin(theta), 0,
		0, 1, 0, 0,
		-sin(theta), 0, cos(theta), 0,
		0, 0, 0, 1;
	Punto = Ry * Punto;
	cop_x = Punto.x();
	cop_y = Punto.y();
	cop_z = Punto.z();
}

void RotacionenX(float theta)
{
	Punto << cop_x, cop_y, cop_z, 1;
	Rx << 1, 0, 0, 0,
		0, cos(theta), -sin(theta), 0,
		0, sin(theta), cos(theta), 0,
		0, 0, 0, 1;
	Punto = Rx * Punto;
	cop_x = Punto.x();
	cop_y = Punto.y();
	cop_z = Punto.z();
	Punto << pref_x, pref_y, pref_z, 1;
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'A':
	case 'a':
		RotacionenY(0.1);
		break;
	case 'D':
	case 'd':
		RotacionenY(-0.1);
		break;
	case 'W':
	case 'w':
		RotacionenX(-0.1);
		break;
	case 'S':
	case 's':
		RotacionenX(0.1);
		break;
	case 'F':
	case 'f':
		cop_x = cop_x - 0.1;
		break;
	case 'H':
	case 'h':
		cop_x = cop_x + 0.1;
		break;
	case 'T': //Zoom +
	case 't':
		cop_z = cop_z - 0.1;
		break;
	case 'G': //Zoom -
	case 'g':
		cop_z = cop_z + 0.1;
		break;
	case 'Z': //Malla
	case 'z':
		M = !M;
		break;
	case 'X': //Ejes
	case 'x':
		E = !E;
		break;
	case 'Q': //Figura
	case 'q':
		Q = !Q;
		break;
	case 'N': //Figura
	case 'n':
		N = !N;
		break;
	}
	glutPostRedisplay();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// glFrustum (-1.0, 1.0, -1.0, 1.0, 1.5, 20.0);
	gluPerspective(60.0, (float)w / (float)h, 1.0, 300.0);
	glFlush();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(xWmax, yWmax);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	init();
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutTimerFunc(1, OnTimerGL, 1);
	glutReshapeFunc(reshape);
	glutMainLoop();
	return 0;
}
