
Room Impulse Response Generator + Mixing

- ��������(20180122)
  1. genNxN.m ����
     �������� mixingNxN.m���� �������� �� RIR�� ������ �ӽ÷� ���� �� �ͽ��ϰ� ���� �������� �ʾ�����
     generateRIR�� ���� RIR ������ �����Ͽ� ������� �����ϵ��� ��(�����ϱ�����)
     

- ��������(20171124)
  1. ���� rir.m(Stephen G. McGovern)�� ����
  2. Original Image method �� ������ �ڵ�(rir_generator)�� ����
     Image method for efficiently simulating small-room acoustics(Jont B. Allen and David A. Berkley, 1978)
  3. ���� ��ĺ��� RIR ������ �ð��� �ſ� ���� �ɸ��Ƿ� �̸� RIR�� .mat ���Ϸ� ���� ��
     �Էµ� room dimension, mic location, azimuth, elevation, RT60�� �´� RIR�� �ҷ��� mixing ����


- mixing_main.m
  1. preparation �κа� mixing �κ��� ����
  2. preparation�� �Էµ� ����ȯ�� ������ ���� �̸� RIR�� �����ϴ� �κ�
  3. 1024kHz sampling���� �����ؾ� �� (mixing �� 1024kHz�� upsampling �� �� RIR�� convolution -> down sampling ����)
  4. ���� �Էµ� configuration�� �´� RIR�� �̹� ����Ǿ� �ִ� ���� �Ѿ��, ���� configuration�� RIR ���� ���� �� ����
  5. mixing �κ��� ������ ���� �����ϸ�, RIR ������ ���ƴٸ� �� �κи� ���� ���ȯ�濡 �°� �����ؼ� ���

- generateRIR.m
  1. Original image method�� �ԷµǴ� ����ȯ�濡 �°� �����ǵ��� ��
  2. distance, azimuth, elevation�� ������ �Է��ϸ� �Է� ������ŭ �����ϰ� �� 

- mixingNxN.m
  1. ������ ���� �����ϳ� RIR ���� ��ĸ� Original image method �� ���·� ������
  2. �� �κи� �ڽ��� ȯ�濡 �°� �Ķ���͸� �Է��ص� ������ ������
     (RIR�� �̸� ����Ǿ� ���� �ʴٸ� ���� ����ϱ� ������ �ð��� �ſ� �����ɸ�)

- generateRIR.m �� genNxN.m ���� ���� rir_generator �Լ�
  1. RIR-Generator-master ���� ���� rir_generator.cpp�� Matlab���� ��� �����ϵ��� MEX bulid�� �Ͽ� �����
  2. ���� RIR ���� ������ rir_generator_original.cpp ������
     (������ ������ copyright ������ ��µǾ� command â �ȿ��� �����Ȳ�� �ľ��� ����� ������
      rir_generator.cpp������ ������ �ּ�ó����)
  3. �ڽ��� PC ȯ�濡�� ���� �߻� �� RIR-Generator-master ���� �ȿ��� rir_generator.cpp�� MEX build�� ���� �����ؼ� ����ؾ� ��
  4. Matlab command�� 'mex rir_generator.cpp' �� �Է��ϸ� mex ������ ������
  5. mex ������ generatorRIR.m�� �ִ� ���� ���� ����(���� mex ������ ����)

2020.07.20
�������� ���峪�� �������� �� ��ġ�� share��.
