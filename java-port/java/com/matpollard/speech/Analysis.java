package com.matpollard.speech;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class Analysis {

  public static void main(String... aArgs) throws IOException{
    Analysis analysis = new Analysis();
    byte[] bytes = analysis.readSmallBinaryFile(FILE_NAME);
    log("Small - size of file read in:" + bytes.length);


    int[] audioData = new int[bytes.length / 4];
    //int nlengthInSamples = bytes.length / 2;
	//ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer().get(shorts);

    //analysis.bytesToShort(bytes, audioData);
    analysis.bytesToInt(bytes, audioData);

    for (int i=0;i<audioData.length;i++) {
    	log(audioData[i]);
    }


    //analysis.writeSmallBinaryFile(bytes, OUTPUT_FILE_NAME);
  }


  public void bytesToShort(byte[] bytes, short[] shorts) {
  	ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer().get(shorts);
  }

  public void bytesToLong(byte[] bytes, long[] longs) {
  	ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asLongBuffer().get(longs);
  }

  public void bytesToInt(byte[] bytes, int[] ints) {
  	ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asIntBuffer().get(ints);
  }

  final static String FILE_NAME = "MES101.EXC";
  final static String OUTPUT_FILE_NAME = "MES101.EXC.OUT";
  
  byte[] readSmallBinaryFile(String aFileName) throws IOException {
    Path path = Paths.get(aFileName);
    return Files.readAllBytes(path);
  }
  
  void writeSmallBinaryFile(byte[] aBytes, String aFileName) throws IOException {
    Path path = Paths.get(aFileName);
    Files.write(path, aBytes); //creates, overwrites
  }
  
  private static void log(Object aMsg){
    System.out.println(String.valueOf(aMsg));
  }
  
}