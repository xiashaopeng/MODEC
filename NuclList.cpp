#include "NuclList.h"

NuclLibrary::NuclLibrary() // construct function for initializing vector 'nuclide_index_vector_'
{
	effective_fission_yields_vector_.resize(30);
	
	nuclide_index_vector_.resize(nuclide_list_[nuclide_number_ - 1] + 1);

	for (int i = 0; i<nuclide_number_; ++i)
	{
		nuclide_index_vector_[nuclide_list_[i]] = i;
	}
	
	nuclide_library_vector_.resize(11);
	for (int i = 0; i < 10; ++i)
	{
		nuclide_library_vector_[i].resize(nuclide_number_); // ��ѡ���ڴ˴��޸���Ҫ�洢�ĺ�����Ϣ
									  /* ��ʱѡ��洢�ĺ�����Ϣ���£�
									  0: N       ����Ũ�ȣ���λ����ν����Ҫͳһ
									  1: lamda   ˥�䳣�������ڷ����Ի�ȼ��� /s
									  2: FSF     �Է��ѱ�ݶ���ڼ����Է��ѱ����ӳ�
									  3: Q       ˥����ϵ��������˥���ȼ���
									  4: AMPC   ORIGEN���ݿ��ж��Ե�λ m3-air/Bq
									  5: WMPC   ORIGEN���ݿ��ж��Ե�λ m3-water/Bq

									  6: tot-cap ���ص��ܵ����ӷ�Ӧ���棬���ڹ��ʵ�ͨ����ת������
									  7: fission ���ص��ѱ䷴Ӧ���棬���ڹ��ʵ�ͨ����ת������
									  8: Radiative capture ���������棬���ڹ��ʵ�ͨ����ת������
									  9: Neutron Product XS ���Ӳ������棬ÿ�κ˷�Ӧ���������Ӹ������������ӵ�ʧƽ�����
									  */
	}
	// nuclide_library_vector_[10] ��ʾ���Բ�������λSv/Bq
	nuclide_library_vector_[10] = { 0.0, 0.0, 1.8e-11, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 2.8e-11,
		0.0, 0.0, 1.1e-09, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0,
		5.8e-10, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 3.2e-09,
		0.0, 4.3e-10, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 2.2e-09, 3.5e-09,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 1.6e-10, 5.6e-10, 0.0,
		2.4e-09, 2.4e-10, 0.0, 0.0,
		0.0, 0.0, 1.3e-10, 0.0, 0.0,
		0.0, 9.3e-10, 0.0, 1.2e-10,
		0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 6.2e-09,
		0.0, 4.3e-10, 2.5e-10, 8.4e-11,
		0.0, 1.9e-10, 0.0, 0.0, 0.0,
		7.1e-10, 0.0, 1.6e-09, 0.0, 0.0,
		3.5e-10, 2.4e-09, 0.0, 0.0,
		1.5e-09, 0.0, 5.4e-10, 1.7e-09,
		8.2e-11, 0.0, 5.8e-09, 1.5e-10,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		2e-09, 1.8e-11, 0.0, 0.0, 0.0,
		0.0, 0.0, 2e-10, 6.1e-11, 0.0,
		3.8e-11, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 1.8e-09, 3e-11, 7.1e-10,
		0.0, 2.5e-10, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 3.3e-10, 0.0, 0.0,
		0.0, 1.8e-09, 1.1e-07, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		1e-09, 2.5e-09, 2.1e-10, 7.4e-10,
		2.4e-11, 0.0, 3.4e-09, 1.7e-12,
		7.4e-11, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 8.6e-10, 8.7e-10,
		0.0, 6.3e-11, 0.0, 0.0, 0.0,
		1.5e-10, 0.0, 1.8e-10, 3e-09,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0,0.0,1.2e-10,0.0,0.0,
		3.4e-10,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,7.9e-11,0.0,
		3.9e-09,0.0,0.0,0.0,3.1e-11,
		3.3e-10,0.0,0.0,2.4e-10,1.4e-09,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,1.2e-09,1.9e-10,
		1e-10,0.0,3.1e-11,0.0,1.1e-09,0.0,
		2.6e-10,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,1e-10,6.5e-11,1.3e-09,
		2.4e-10,0.0,1.2e-11,0.0,0.0,0.0,
		0.0,0.0,4.6e-11,0.0,0.0,3.3e-10,
		0.0,1.2e-10,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		5.7e-11,4.6e-10,1.8e-09,2.6e-10,
		1.3e-09,0.0,0.0,1.6e-09,4e-10,2.1e-10,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,2.1e-10,2.8e-11,
		0.0,2.6e-09,0.0,
		0.0,
		0.0,
		0.0,
		2.9e-09,
		0.0,
		0.0,
		2.7e-11,
		5.3e-11,
		0.0,
		4.7e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		7.9e-11,
		4.6e-10,
		9.6e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		3.1e-11,
		1.1e-10,
		0.0,
		5.4e-10,
		0.0,
		4.3e-11,
		8.8e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		5e-11,
		5.4e-11,
		0.0,
		1.9e-09,
		2.8e-09,
		0.0,
		2.8e-09,
		0.0,
		1.5e-09,
		9e-11,
		4.7e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		6.1e-09,
		4.9e-10,
		0.0,
		5.6e-10,
		6.1e-12,
		0.0,
		0.0,
		3e-11,
		0.0,
		2.6e-09,
		2.8e-08,
		6.5e-10,
		4.3e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		9.6e-10,
		5.5e-10,
		0.0,
		1.3e-09,
		0.0,
		0.0,
		2.7e-09,
		1.7e-10,
		2.4e-09,
		1.1e-11,
		4.9e-10,
		1.2e-09,
		0.0,
		8.1e-11,
		4.6e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		8.6e-10,
		0.0,
		4.5e-10,
		7.9e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		1.1e-09,
		0.0,
		9.5e-10,
		0.0,
		2.1e-09,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.7e-10,
		1.2e-09,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		1.2e-10,
		1.7e-09,
		0.0,
		5.8e-10,
		5.6e-10,
		1.1e-09,
		6.8e-11,
		0.0,
		1.1e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.2e-10,
		0.0,
		0.0,
		3.1e-09,
		1.1e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		6e-10,
		0.0,
		4.1e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		5.5e-11,
		1.8e-10,
		5.6e-10,
		1.1e-09,
		6.8e-11,
		5.5e-10,
		2e-09,
		6.4e-10,
		2.2e-11,
		0.0,
		1.9e-11,
		0.0,
		0.0,
		0.0,
		8e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		1.5e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		7.3e-10,
		0.0,
		2.6e-10,
		7e-09,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		5.1e-10,
		6.6e-11,
		7.1e-10,
		5.5e-10,
		2.2e-10,
		2.6e-09,
		1.2e-09,
		0.0,
		3.8e-12,
		0.0,
		0.0,
		3.7e-10,
		0.0,
		0.0,
		1.6e-10,
		2.4e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		9.4e-10,
		9.4e-11,
		0.0,
		1.9e-10,
		0.0,
		0.0,
		0.0,
		3.7e-11,
		0.0,
		0.0,
		5.5e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		4.3e-11,
		4.7e-10,
		0.0,
		3.2e-11,
		1.5e-09,
		0.0,
		0.0,
		0.0,
		2.3e-09,
		0.0,
		0.0,
		0.0,
		2.8e-09,
		1.3e-09,
		0.0,
		4.3e-10,
		0.0,
		0.0,
		0.0,
		6e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		6.2e-11,
		0.0,
		2e-09,
		0.0,
		0.0,
		0.0,
		0.0,
		2.5e-08,
		2.3e-08,
		0.0,
		1.4e-09,
		3.3e-09,
		0.0,
		2.8e-10,
		2.8e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		6.6e-11,
		2.9e-10,
		0.0,
		1e-11,
		0.0,
		0.0,
		2.8e-11,
		0.0,
		4.1e-09,
		3.2e-08,
		8.6e-11,
		0.0,
		6.4e-11,
		3.1e-11,
		1.2e-10,
		0.0,
		0.0,
		0.0,
		4.7e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.3e-11,
		0.0,
		7.3e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		7.1e-10,
		0.0,
		0.0,
		3.4e-10,
		0.0,
		2.3e-10,
		3.8e-10,
		0.0,
		2.1e-09,
		3.8e-11,
		0.0,
		3.1e-09,
		0.0,
		4.7e-09,
		2e-10,
		0.0,
		1.5e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.4e-11,
		1.8e-11,
		0.0,
		2.1e-10,
		8e-11,
		1.4e-11,
		1.2e-09,
		0.0,
		1.7e-09,
		0.0,
		0.0,
		2.5e-09,
		0.0,
		1.1e-09,
		2.4e-09,
		3.6e-11,
		1.7e-09,
		7.6e-10,
		3.3e-11,
		4.2e-10,
		0.0,
		9.1e-11,
		0.0,
		1e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		4.3e-10,
		2.3e-09,
		0.0,
		4.4e-09,
		1.4e-09,
		0.0,
		0.0,
		8.7e-10,
		0.0,
		1.7e-10,
		2.3e-09,
		0.0,
		6.3e-11,
		3e-09,
		0.0,
		8.7e-11,
		1.9e-09,
		3.8e-09,
		7.2e-11,
		2.8e-10,
		1.1e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		8.2e-11,
		0.0,
		2.1e-10,
		1.3e-08,
		1.5e-08,
		2.9e-08,
		0.0,
		4.6e-11,
		1.1e-07,
		2e-09,
		0.0,
		2.2e-08,
		2.9e-10,
		2.2e-10,
		4.3e-09,
		0.0,
		1.1e-10,
		0.0,
		9.3e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.4e-11,
		0.0,
		6e-11,
		2.8e-11,
		5.8e-11,
		5e-10,
		0.0,
		1.9e-08,
		2e-11,
		2e-09,
		1.9e-11,
		3e-09,
		0.0,
		1.3e-08,
		9.2e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.7e-09,
		0.0,
		0.0,
		4.5e-10,
		4.9e-12,
		0.0,
		1.5e-09,
		5.4e-10,
		0.0,
		0.0,
		4.3e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		1.2e-10,
		2.6e-09,
		7e-11,
		3.5e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		3e-11,
		0.0,
		8.1e-11,
		1.1e-09,
		0.0,
		2e-09,
		3.6e-10,
		1.8e-10,
		5.6e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.5e-09,
		7.9e-10,
		0.0,
		2.5e-11,
		5.4e-10,
		0.0,
		2.6e-10,
		0.0,
		0.0,
		7.1e-10,
		0.0,
		1.1e-09,
		5.2e-09,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		3.1e-11,
		0.0,
		0.0,
		1.3e-09,
		1.7e-11,
		1.2e-09,
		5e-11,
		0.0,
		3.9e-10,
		0.0,
		3.3e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		8.3e-12,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		1.1e-09,
		0.0,
		1.2e-10,
		0.0,
		3e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		3.6e-11,
		2.3e-10,
		9.7e-10,
		1.1e-10,
		9e-10,
		2.6e-10,
		2.7e-09,
		1.7e-09,
		9.9e-10,
		2.6e-10,
		7.3e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		2.1e-10,
		5.4e-08,
		4.9e-08,
		0.0,
		0.0,
		0.0,
		9.8e-11,
		0.0,
		7.4e-10,
		0.0,
		2.9e-11,
		2.5e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		7.5e-10,
		1.3e-09,
		4.4e-10,
		1.3e-09,
		1e-10,
		1.3e-09,
		3.8e-10,
		0.0,
		1.4e-09,
		5e-10,
		0.0,
		2e-09,
		0.0,
		3.2e-10,
		2.2e-09,
		6e-10,
		9.4e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		9.6e-10,
		6.1e-10,
		5.6e-08,
		4.5e-10,
		0.0,
		2e-10,
		4.1e-08,
		2.7e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		4.9e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		3.4e-10,
		0.0,
		2.5e-10,
		6.5e-10,
		0.0,
		2.1e-10,
		1.2e-09,
		1.7e-10,
		3.4e-11,
		1.1e-09,
		0.0,
		0.0,
		1.6e-09,
		7.2e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		1.3e-10,
		0.0,
		6.1e-11,
		0.0,
		1e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		1.1e-10,
		0.0,
		1.6e-09,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		7.9e-12,
		0.0,
		0.0,
		0.0,
		1.3e-11,
		0.0,
		3.3e-12,
		2.6e-11,
		0.0,
		0.0,
		9.5e-12,
		1.6e-11,
		0.0,
		1.4e-09,
		2e-09,
		8.3e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		8e-11,
		0.0,
		0.0,
		0.0,
		1.9e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		3.7e-10,
		0.0,
		3.6e-10,
		1e-09,
		0.0,
		2.8e-10,
		5.6e-10,
		0.0,
		0.0,
		1.3e-09,
		0.0,
		1.1e-10,
		1.7e-09,
		3.1e-10,
		9.5e-10,
		6.7e-12,
		0.0,
		7.1e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		4.4e-10,
		0.0,
		0.0,
		8.8e-11,
		4.6e-10,
		0.0,
		9.9e-10,
		6.7e-10,
		0.0,
		1.3e-09,
		0.0,
		2.6e-10,
		2.7e-10,
		5.3e-10,
		0.0,
		1.8e-09,
		1.7e-10,
		5.3e-10,
		1.7e-09,
		4.8e-10,
		0.0,
		1e-09,
		2.3e-10,
		0.0,
		4.1e-10,
		0.0,
		0.0,
		8.1e-11,
		0.0,
		4.7e-09,
		0.0,
		1.2e-09,
		0.0,
		1.7e-10,
		1.1e-09,
		3e-09,
		1.1e-10,
		0.0,
		6.5e-11,
		8.4e-10,
		5.4e-11,
		0.0,
		1.5e-09,
		1.2e-11,
		1.3e-09,
		2.2e-10,
		0.0,
		7.6e-11,
		0.0,
		0.0,
		0.0,
		0.0,
		4.4e-10,
		0.0,
		0.0,
		6.3e-10,
		2.1e-09,
		0.0,
		4.2e-10,
		1.4e-09,
		2.7e-10,
		0.0,
		1e-09,
		1.5e-09,
		0.0,
		1.5e-09,
		2.2e-09,
		5.1e-12,
		1.4e-09,
		3e-11,
		7.8e-10,
		5.6e-10,
		0.0,
		0.0,
		5.1e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		1.8e-11,
		0.0,
		0.0,
		5.7e-10,
		9.6e-11,
		0.0,
		8.1e-10,
		2.4e-09,
		2.6e-10,
		4.9e-10,
		6.3e-10,
		2.4e-10,
		0.0,
		1.2e-09,
		0.0,
		0.0,
		1.4e-09,
		3.1e-10,
		0.0,
		2.7e-10,
		1.3e-09,
		2.1e-09,
		0.0,
		0.0,
		7.6e-10,
		1.2e-10,
		0.0,
		3.4e-10,
		0.0,
		3.1e-11,
		4.5e-10,
		0.0,
		0.0,
		6.3e-10,
		0.0,
		4e-10,
		8.4e-11,
		0.0,
		3.9e-11,
		0.0,
		1.2e-09,
		1.3e-10,
		4.2e-10,
		2.5e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		1e-09,
		1.3e-09,
		4.4e-10,
		6.8e-11,
		1.1e-09,
		8.2e-11,
		4e-10,
		1.4e-09,
		9.7e-11,
		5.6e-10,
		0.0,
		2.3e-10,
		4.7e-10,
		0.0,
		0.0,
		3.1e-11,
		0.0,
		0.0,
		0.0,
		5.4e-10,
		0.0,
		0.0,
		0.0,
		2e-10,
		9.5e-11,
		4.5e-10,
		0.0,
		1.2e-09,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		4e-10,
		8.8e-09,
		2.4e-10,
		0.0,
		2.8e-10,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		5.7e-11,
		6.9e-07,
		1.8e-10,
		6e-09,
		1.4e-10,
		9e-10,
		1.9e-09,
		1.3e-09,
		0.0,
		0.0,
		1.3e-09,
		1.5e-08,
		0.0,
		2.6e-10,
		0.0,
		2e-10,
		1.1e-10,
		0.0,
		1.1e-10,
		0.0,
		0.0,
		1.2e-06,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		7.2e-10,
		2.4e-09,
		0.0,
		0.0,
		1e-07,
		6.5e-08,
		9.9e-08,
		2.8e-07,
		8.1e-11,
		6.9e-07,
		7e-10,
		2.4e-08,
		1e-08,
		1.1e-06,
		4.3e-10,
		3.5e-10,
		8.8e-09,
		7.2e-08,
		4.9e-07,
		2.1e-07,
		3.4e-10,
		2.3e-07,
		0.0,
		3.4e-09,
		7.8e-10,
		0.0,
		9.2e-10,
		7.1e-07,
		7.2e-10,
		8.7e-10,
		5.1e-10,
		0.0,
		0.0,
		5.6e-08,
		2.8e-10,
		3.3e-07,
		5.1e-08,
		4.9e-08,
		4.7e-08,
		0.0,
		4.7e-08,
		7.6e-10,
		4.5e-08,
		2.7e-11,
		1.1e-09,
		0.0,
		8.1e-10,
		5.3e-11,
		1.7e-08,
		1.9e-10,
		1.1e-07,
		9.1e-10,
		8e-10,
		8.2e-11,
		0.0,
		0.0,
		8.7e-08,
		1e-10,
		0.0,
		2.3e-07,
		2.5e-07,
		2.5e-07,
		4.8e-09,
		2.4e-07,
		8.5e-11,
		2.4e-07,
		7.2e-10,
		3.3e-09,
		0.0,
		2.4e-10,
		5.8e-10,
		2e-07,
		3e-10,
		1.9e-07,
		2e-07,
		4.6e-10,
		2.9e-11,
		6.2e-11,
		5.8e-11,
		3.4e-11,
		0.0,
		7.6e-09,
		9.1e-10,
		1.2e-08,
		1.5e-07,
		1.2e-07,
		2.1e-07,
		2.1e-07,
		1.9e-07,
		7.7e-07,
		3.1e-11,
		4.4e-06,
		0.0,
		5.7e-10,
		4.8e-10,
		3.5e-07,
		0.0,
		0.0,
		9.7e-10,
		1.4e-10,
		0.0,
		3.3e-09,
		2.8e-08,
		3.5e-07,
		1.6e-07,
		3.6e-07,
		9e-08,
		1.4e-09,
		4e-07,
		0.0,
		1.7e-10,
		0.0,
		6.1e-09,
		2.8e-08,
		4.2e-09,
		0.0
	};
};

vector<int> NuclLibrary::GetEleIndex(int Z) // �ṩԭ����������������ͬλ��
{
	int ID = Z * 10000; // ��С��ID��������ӦΪ1
	vector<int> IDList;
	int ListEnd;

	for (int i = 0; i < nuclide_number_; ++i)
	{
		if (nuclide_list_[i] >= ID && nuclide_list_[i] < ID + 10000)
		{
			IDList.push_back(i);
		}
		if (nuclide_list_[i] >= ID + 10000)
		{
			ListEnd = i - 1;
			break;
		}
	}
	return IDList;
}