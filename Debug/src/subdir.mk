################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Body.cpp \
../src/ElasticBody.cpp \
../src/IntegratorCompareVisitor.cpp \
../src/IntegratorGNUPlotVisitor.cpp \
../src/MBSystem.cpp \
../src/ODEOrder2.cpp 

OBJS += \
./src/Body.o \
./src/ElasticBody.o \
./src/IntegratorCompareVisitor.o \
./src/IntegratorGNUPlotVisitor.o \
./src/MBSystem.o \
./src/ODEOrder2.o 

CPP_DEPS += \
./src/Body.d \
./src/ElasticBody.d \
./src/IntegratorCompareVisitor.d \
./src/IntegratorGNUPlotVisitor.d \
./src/MBSystem.d \
./src/ODEOrder2.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -I/usr/local/include/eigen3 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


