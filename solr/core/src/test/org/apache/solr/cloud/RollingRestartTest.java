package org.apache.solr.cloud;

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import org.apache.solr.common.cloud.SolrZkClient;
import org.apache.solr.common.params.CollectionParams;
import org.apache.zookeeper.KeeperException;
import org.junit.After;
import org.junit.Before;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

public class RollingRestartTest extends AbstractFullDistribZkTestBase {
  public static Logger log = LoggerFactory.getLogger(ChaosMonkeyNothingIsSafeTest.class);

  public RollingRestartTest() {
    fixShardCount = true;
    sliceCount = 2;
    shardCount = 16;
  }

  @Before
  @Override
  public void setUp() throws Exception {
    super.setUp();
    System.setProperty("numShards", Integer.toString(sliceCount));
    useFactory("solr.StandardDirectoryFactory");
  }

  @Override
  @After
  public void tearDown() throws Exception {
    System.clearProperty("numShards");
    super.tearDown();
    resetExceptionIgnores();
  }

  @Override
  public void doTest() throws Exception {
    waitForRecoveriesToFinish(false);

    restartWithRolesTest();

    waitForRecoveriesToFinish(false);
  }


  public void restartWithRolesTest() throws Exception {
    String leader = OverseerCollectionProcessor.getLeaderNode(cloudClient.getZkStateReader().getZkClient());
    assertNotNull(leader);
    log.info("Current overseer leader = {}", leader);

    cloudClient.getZkStateReader().getZkClient().printLayoutToStdOut();

    int numOverseers = 3;
    List<String> designates = new ArrayList<>();
    List<CloudJettyRunner> overseerDesignates = new ArrayList<>();
    for (int i = 0; i < numOverseers; i++) {
      int n = random().nextInt(shardCount);
      String nodeName = cloudJettys.get(n).nodeName;
      log.info("Chose {} as overseer designate", nodeName);
      invokeCollectionApi(CollectionParams.ACTION, CollectionParams.CollectionAction.ADDROLE.toLower(), "role", "overseer", "node", nodeName);
      designates.add(nodeName);
      overseerDesignates.add(cloudJettys.get(n));
    }

    waitUntilOverseerDesignateIsLeader(cloudClient.getZkStateReader().getZkClient(), designates, 60);

    cloudClient.getZkStateReader().getZkClient().printLayoutToStdOut();

    int numRestarts = 4; // 1 + random().nextInt(5);
    for (int i = 0; i < numRestarts; i++) {
      log.info("Rolling restart #{}", i + 1);
      for (CloudJettyRunner cloudJetty : overseerDesignates) {
        log.info("Restarting {}", cloudJetty);
        chaosMonkey.stopJetty(cloudJetty);
        boolean success = waitUntilOverseerDesignateIsLeader(cloudClient.getZkStateReader().getZkClient(), designates, 60);
        if (!success) {
          leader = OverseerCollectionProcessor.getLeaderNode(cloudClient.getZkStateReader().getZkClient());
          if(leader == null) log.error("NOOVERSEER election queue is :"+ OverseerCollectionProcessor.getSortedElectionNodes(cloudClient.getZkStateReader().getZkClient()));
          fail("No overseer designate as leader found after restart #" + (i + 1) + ": " + leader);
        }
        assertTrue("Unable to restart (#"+i+"): " + cloudJetty, 
                   chaosMonkey.start(cloudJetty.jetty));
        success = waitUntilOverseerDesignateIsLeader(cloudClient.getZkStateReader().getZkClient(), designates, 60);
        if (!success) {
          leader = OverseerCollectionProcessor.getLeaderNode(cloudClient.getZkStateReader().getZkClient());
          if(leader == null) log.error("NOOVERSEER election queue is :"+ OverseerCollectionProcessor.getSortedElectionNodes(cloudClient.getZkStateReader().getZkClient()));
          fail("No overseer leader found after restart #" + (i + 1) + ": " + leader);
        }
      }
    }

    leader = OverseerCollectionProcessor.getLeaderNode(cloudClient.getZkStateReader().getZkClient());
    assertNotNull(leader);
    log.info("Current overseer leader (after restart) = {}", leader);

    cloudClient.getZkStateReader().getZkClient().printLayoutToStdOut();
  }

  static boolean waitUntilOverseerDesignateIsLeader(SolrZkClient testZkClient, List<String> overseerDesignates, int timeoutInSeconds) throws KeeperException, InterruptedException {
    long now = System.nanoTime();
    long timeout = now + TimeUnit.NANOSECONDS.convert(timeoutInSeconds, TimeUnit.SECONDS);
    boolean firstTime = true;
    int stableCheckTimeout = 2000;
    while (System.nanoTime() < timeout) {
      String newLeader = OverseerCollectionProcessor.getLeaderNode(testZkClient);
      if (!overseerDesignates.contains(newLeader)) {
        Thread.sleep(500);
      } else {
        if (firstTime)  {
          firstTime = false;
          Thread.sleep(stableCheckTimeout);
        } else  {
          return true;
        }
      }
    }
    return false;
  }
}
